#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <stdio.h>
#include <map>

//#include "moab/CartVect.hpp"
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/ElemEvaluator.hpp"
#include "moab/CN.hpp"
#include "moab/SpatialLocator.hpp"
#include "moab/Util.hpp"
#include "moab/GeomUtil.hpp"
#include "MBTagConventions.hpp"

moab::Core mbi;
moab::Tag result_tag;
moab::Tag error_tag;
moab::Tag id_tag;
moab::Tag sdr_tag;
moab::Tag sdr_err_tag;
moab::Tag p_src_err_tag;
moab::Tag p_src_tag;
moab::Tag ap_flux_tag;
moab::Tag sq_p_src_err_tag2;

struct Tet_info{
       moab::EntityHandle eh;
       std::vector<double> result;
       double vol;
};

struct CartVect{
    double coords[3];
    inline CartVect() {}

    inline CartVect( double tx, double ty, double tz ) { set(tx,ty,tz); }

    inline CartVect( const CartVect& other ) { set( other.coords); }

    inline void set( double tx, double ty, double tz )
      { coords[0] = tx; coords[1] = ty; coords[2] = tz; }

    inline void set( const double* c )
      { coords[0] = c[0]; coords[1] = c[1]; coords[2] = c[2]; }

    inline double x() const { return coords[0]; }
    inline double y() const { return coords[1]; }
    inline double z() const { return coords[2]; }

    inline CartVect& operator+=( const CartVect& other )
    {
      coords[0] += other.coords[0];
      coords[1] += other.coords[1];
      coords[2] += other.coords[2];
      return *this;
    }

    inline CartVect& operator-=( const CartVect& other )
    {
      coords[0] -= other.coords[0];
      coords[1] -= other.coords[1];
      coords[2] -= other.coords[2];
      return *this;
    }

    inline CartVect& operator*=( const CartVect& other );

    inline double lensqr() const;

    inline double len() const;

    inline CartVect operator~( ) const;


    inline CartVect& operator*=( double a )
    {
      coords[0] *= a;
      coords[1] *= a;
      coords[2] *= a;
      return *this;
    }

    inline CartVect& operator/=( double a )
    {
      coords[0] /= a;
      coords[1] /= a;
      coords[2] /= a;
      return *this;
    }
};

inline CartVect operator-( const CartVect& v1, const CartVect& v2 )
{
  CartVect rval(v1);
  rval -= v2;
  return rval;
}

inline double operator%( const CartVect& v1, const CartVect& v2 )
{
  return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline CartVect operator*( const CartVect& v1, const CartVect& v2 )
{
  return CartVect( v1.y() * v2.z() - v1.z() * v2.y(),
                   v1.z() * v2.x() - v1.x() * v2.z(),
                   v1.x() * v2.y() - v1.y() * v2.x() );
}

static double tet_volume( const CartVect& v0,
                                 const CartVect& v1,
                                 const CartVect& v2,
                                 const CartVect& v3 )
{
  return 1./6. * ( ((v1 - v0) * (v2 - v0)) % (v3 - v0) );
}

int num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mbi.tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}


moab::ErrorCode get_mesh_elements(std::string filename,
                                  std::string result_tag_name,
                                  std::map<int, Tet_info> &tet_map,
                                  moab::EntityHandle &fileset,
                                  int &num_e_groups){

moab::ErrorCode rval;
moab::Range ves;
moab::Range::iterator it;

// Load mesh from file into fileset
rval = mbi.create_meshset(moab::MESHSET_SET, fileset); MB_CHK_ERR(rval);
MB_CHK_SET_ERR(rval, "Error creating meshset.");
rval = mbi.load_file(filename.c_str(), &fileset);
MB_CHK_SET_ERR(rval, "Error loading file.");

// Get ID tag
//rval = mbi.tag_get_handle( GLOBAL_ID_TAG_NAME,
rval = mbi.tag_get_handle( "GLOBAL_ID",
                           1, 
                           moab::MB_TYPE_INTEGER,
                           id_tag,
                           moab::MB_TAG_DENSE );
MB_CHK_SET_ERR(rval, "Error getting ID handle.");
// Get result tag
std::cout << "result tag name: " << result_tag_name << std::endl;
rval = mbi.tag_get_handle(result_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           result_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//                           moab::MB_TAG_DENSE);//|moab::MB_TAG_CREAT);
//MB_CHK_SET_ERR(rval, "Error getting result tag handle.");
num_e_groups = num_groups(result_tag);
std::cout << "num e groups " << num_e_groups << std::endl;
std::vector<double> result(num_e_groups, 0);

// Get all 3D elements in fileset 
ves.clear();
//rval = mbi.get_entities_by_dimension(fileset, 3, ves);
rval = mbi.get_entities_by_type(fileset, moab::MBTET, ves);
MB_CHK_SET_ERR(rval, "Error getting tets.");
std::cout << "num ves " << ves.size() << std::endl;
int tet_id;
for (it = ves.begin(); it != ves.end(); ++it){
  // get the id tag on the ve
  rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id );
  MB_CHK_SET_ERR(rval, "Error getting id tag.");

  // get the result tag on the ve 
  rval = mbi.tag_get_data(result_tag, &(*it), 1, &(result[0]));
  MB_CHK_SET_ERR(rval, "Error getting result tag.");

  // get the volume of the ve
  moab::Range vertices;
  double coord[3];
  rval = mbi.get_connectivity(&(*it), 1, vertices); 
  MB_CHK_SET_ERR(rval, "Error getting vertices.");
  // for each vertex, get coords
  CartVect coords[4];
  int i=0;
  for (moab::Range::iterator it = vertices.begin(); it != vertices.end(); it++){
    rval = mbi.get_coords(&(*it), 1 , &(coord[0]));
    MB_CHK_SET_ERR(rval, "Error getting coords.");
    coords[i].coords[0] = coord[0];
    coords[i].coords[1] = coord[1];
    coords[i].coords[2] = coord[2];
    ++i;
   // std::cout << "coords " << coord[0] << " " << coord[1] << " " << coord[2] << std::endl;
  }
  double ve_vol = fabs(tet_volume( coords[0], coords[1], coords[2], coords[3] ));
  //double ve_vol = tet_volume( coords[0], coords[1], coords[2], coords[3] );
  //std::cout << "tet vol" << tet_id << " " << ve_vol << std::endl;

  // Add entity handle, result, and vol to tet id map 
  tet_map[tet_id].eh = *it;
  tet_map[tet_id].result = result;
  tet_map[tet_id].vol = ve_vol;
}

return moab::MB_SUCCESS;
}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::EntityHandle fileset;

std::map<int, Tet_info> adj_flux_map;
std::map<int, Tet_info> p_src_map;
std::map<int, Tet_info> sq_err_p_src_map;

// Get p src err tag
std::string p_src_err_tag_name ("p_src_err");
rval = mbi.tag_get_handle(p_src_err_tag_name.c_str(),
                           //moab::MB_TAG_VARLEN,
                           42,
                           moab::MB_TYPE_DOUBLE,
                           p_src_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string sq_p_src_err_tag_name ("sq_p_src_err2");
rval = mbi.tag_get_handle(sq_p_src_err_tag_name.c_str(),
                           //moab::MB_TAG_VARLEN,
                           42,
                           moab::MB_TYPE_DOUBLE,
                           sq_p_src_err_tag2,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string ap_flux_tag_name ("ap_flux");
rval = mbi.tag_get_handle(ap_flux_tag_name.c_str(),
                           //moab::MB_TAG_VARLEN,
                           217,
                           moab::MB_TYPE_DOUBLE,
                           ap_flux_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
// Get sdr tag
std::string sdr_tag_name ("sdr");
rval = mbi.tag_get_handle(sdr_tag_name.c_str(),
                           //moab::MB_TAG_VARLEN,
                           42,
                           moab::MB_TYPE_DOUBLE,
                           sdr_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
// Get sdr err tag
std::string sdr_err_tag_name ("sdr_err");
rval = mbi.tag_get_handle(sdr_err_tag_name.c_str(),
                           //moab::MB_TAG_VARLEN,
                           42,
                           moab::MB_TYPE_DOUBLE,
                           sdr_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
// Get p src tag
std::string p_src_tag_name2 ("p_src");
rval = mbi.tag_get_handle(p_src_tag_name2.c_str(),
                           //moab::MB_TAG_VARLEN,
                           42,
                           moab::MB_TYPE_DOUBLE,
                           p_src_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);

// get adjoint photon flux results
//std::cout << "Adjoint photon flux file: " << argv[1] << std::endl;
adj_flux_map.clear();
std::string adj_flux_tag_name = "flux";
int num_e_groups_flux;
rval = get_mesh_elements(argv[1], adj_flux_tag_name, adj_flux_map, fileset, num_e_groups_flux);
MB_CHK_SET_ERR(rval, "Error getting adj flux mesh file");

// get forward photon source results
//std::cout << "Photon source file: " << argv[1] << std::endl;
p_src_map.clear();
//std::string p_src_tag_name = "source_density";
std::string p_src_tag_name = "est_p_src";
int num_e_groups_src;
rval = get_mesh_elements(argv[1], p_src_tag_name, p_src_map, fileset, num_e_groups_src);
MB_CHK_SET_ERR(rval, "Error getting source mesh file");

// get photon source error
//std::cout << "Photon source error file: " << argv[1] << std::endl;
sq_err_p_src_map.clear();
std::string sq_err_p_src_tag_name = "sq_err_p_src";
int num_e_groups_p_err;
rval = get_mesh_elements(argv[1], sq_err_p_src_tag_name, sq_err_p_src_map, fileset, num_e_groups_p_err);
MB_CHK_SET_ERR(rval, "Error getting sq err p src mesh file");

// initialize sdr, sdr err, p src err maps
//map<X, Y> mp2(mp1);
std::map<int, Tet_info> sq_err_sdr_map =p_src_map;
std::map<int, Tet_info> rel_err_sdr_map = p_src_map;
std::map<int, Tet_info> rel_err_p_src_map = p_src_map;
std::map<int, Tet_info> sdr_contributions_map = p_src_map;

// create final meshset
moab::EntityHandle result_meshset;
rval = mbi.create_meshset(moab::MESHSET_SET, result_meshset); MB_CHK_ERR(rval);
MB_CHK_SET_ERR(rval, "Error creating meshset.");


double sdr_contributions = 0.0;
double sdr_at_detector = 0.0;
double sq_err_sdr = 0.0;
moab::EntityHandle tet_id;
moab::EntityHandle p_src_eh;
double vol;
double adj_flux_h;
double p_src_h;
double sq_err_p_src_h;
// for each tet ID, keep running total of the sdr and error
std::map<int, Tet_info>::iterator mit;
std::cout << p_src_map.size() << std::endl;

for(mit = p_src_map.begin(); mit!=p_src_map.end(); ++mit){

  tet_id = mit->first;
  vol = p_src_map[tet_id].vol;
  p_src_eh = p_src_map[tet_id].eh;

  std::vector<double> sdr_result(num_e_groups_src, 0);
  std::vector<double> sq_err_sdr_result(num_e_groups_src, 0);
  std::vector<double> rel_err_sdr_result(num_e_groups_src, 0);
  std::vector<double> rel_err_p_src_result(num_e_groups_src, 0);

  for(int h=0; h <= num_e_groups_src-1; h++){
    adj_flux_h = adj_flux_map[tet_id].result[h];
    p_src_h = p_src_map[tet_id].result[h];
    sq_err_p_src_h = sq_err_p_src_map[tet_id].result[h];

    // relative error in photon source per energy group, per tet
    if (p_src_h == 0){
      rel_err_p_src_result[h] = 0.0;
    }
    else{
      rel_err_p_src_result[h] = (sqrt(sq_err_p_src_h))/p_src_h;
    }
    std::cout << "p src, sq err, rel err " << p_src_h << " " << sq_err_p_src_h << " " << rel_err_p_src_result[h] << std::endl;
    // sdr contribution in each tet, per energy group
    sdr_result[h] = adj_flux_h*p_src_h*vol;

    // accumulate sdr contributions over all energy group, all  tets
    sdr_contributions += sdr_result[h];


    // relative error in sdr contribution per energy group, per tet
    if (sdr_result[h] == 0){
      sq_err_sdr_result[h] = 0.0;
      rel_err_sdr_result[h] = 0.0;
    }
    else{
      sq_err_sdr_result[h] = pow((adj_flux_h*sqrt(sq_err_p_src_h)*vol), 2);
      rel_err_sdr_result[h] = (sqrt(sq_err_sdr_result[h]))/sdr_result[h];
    }
    std::cout << "h sdr, sq err, rel err " << h << " " << sdr_result[h] << " " << sq_err_sdr_result[h] << " " << rel_err_sdr_result[h] << std::endl;
    // accumulate sq. error in sdr over all energy groups, all tets
    //sq_err_sdr += adj_flux_h*sq_err_p_src_h*vol;
    sq_err_sdr += sq_err_sdr_result[h];
    //std::cout << "accum sdr sq err " << sq_err_sdr << std::endl;
  }
  //set the tag vals 
  rval = mbi.tag_set_data(p_src_err_tag, &(p_src_eh), 1, &(rel_err_p_src_result[0]));
  MB_CHK_SET_ERR(rval, "Error setting err in sdr contribution tag val.");

  rval = mbi.tag_set_data(sdr_tag, &(p_src_eh), 1, &(sdr_result[0]));
  MB_CHK_SET_ERR(rval, "Error setting sdr contribution tag val.");

  rval = mbi.tag_set_data(sdr_err_tag, &(p_src_eh), 1, &(rel_err_sdr_result[0]));
  MB_CHK_SET_ERR(rval, "Error setting err in sdr contribution tag val.");


  //add tet to meshset
  rval = mbi.add_entities(result_meshset, &(p_src_eh), 1);
  MB_CHK_SET_ERR(rval, "Error adding EH to result meshset.");
}

//rval = mbi.write_mesh("sdr_err_mesh.h5m");
moab::EntityHandle output_list[] = {result_meshset};
rval = mbi.write_mesh("sdr_err_mesh.h5m", output_list, 1);
MB_CHK_SET_ERR(rval, "Error writing mesh.");

sdr_at_detector = sdr_contributions;
double sdr_err = sqrt(sq_err_sdr);
double sdr_rel_err = (sqrt(sq_err_sdr))/sdr_at_detector;

std::string ts = argv[2];
std::cout << ts << " " << sdr_at_detector << " " << sdr_err << " " << sdr_rel_err << std::endl;
  
}
