#include <iostream>
#include <fstream>
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
moab::Tag total_sdr_tag;
moab::Tag sdr_err_tag;
moab::Tag total_sdr_err_tag;
moab::Tag total_sdr_abs_err_tag;
moab::Tag total_sdr_sq_err_tag;
moab::Tag fom_tag;
moab::Tag total_fom_tag;
moab::Tag p_src_err_tag;
moab::Tag total_p_src_err_tag;

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
//std::cout << "result tag name: " << result_tag_name << std::endl;
rval = mbi.tag_get_handle(result_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           result_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//                           moab::MB_TAG_DENSE);//|moab::MB_TAG_CREAT);
//MB_CHK_SET_ERR(rval, "Error getting result tag handle.");
num_e_groups = num_groups(result_tag);
//std::cout << "num e groups " << num_e_groups << std::endl;
std::vector<double> result(num_e_groups, 0);

// Get all 3D elements in fileset 
ves.clear();
//rval = mbi.get_entities_by_dimension(fileset, 3, ves);
rval = mbi.get_entities_by_type(fileset, moab::MBTET, ves);
MB_CHK_SET_ERR(rval, "Error getting tets.");
//std::cout << "num ves " << ves.size() << std::endl;
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
  }
  double ve_vol = fabs(tet_volume( coords[0], coords[1], coords[2], coords[3] ));

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
                           42,
                           moab::MB_TYPE_DOUBLE,
                           p_src_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string total_p_src_err_tag_name ("total_p_src_err");
rval = mbi.tag_get_handle(total_p_src_err_tag_name.c_str(),
                           1,
                           moab::MB_TYPE_DOUBLE,
                           total_p_src_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
// Get sdr tag
std::string sdr_tag_name ("sdr");
rval = mbi.tag_get_handle(sdr_tag_name.c_str(),
                           42,
                           moab::MB_TYPE_DOUBLE,
                           sdr_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string total_sdr_tag_name ("total_sdr");
rval = mbi.tag_get_handle(total_sdr_tag_name.c_str(),
                           1,
                           moab::MB_TYPE_DOUBLE,
                           total_sdr_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
// Get sdr err tag
std::string sdr_err_tag_name ("sdr_err");
rval = mbi.tag_get_handle(sdr_err_tag_name.c_str(),
                           42,
                           moab::MB_TYPE_DOUBLE,
                           sdr_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string total_sdr_err_tag_name ("total_sdr_err");
rval = mbi.tag_get_handle(total_sdr_err_tag_name.c_str(),
                           1,
                           moab::MB_TYPE_DOUBLE,
                           total_sdr_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string total_sdr_abs_err_tag_name ("total_sdr_abs_err");
rval = mbi.tag_get_handle(total_sdr_abs_err_tag_name.c_str(),
                           1,
                           moab::MB_TYPE_DOUBLE,
                           total_sdr_abs_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string total_sdr_sq_err_tag_name ("total_sdr_sq_err");
rval = mbi.tag_get_handle(total_sdr_sq_err_tag_name.c_str(),
                           1,
                           moab::MB_TYPE_DOUBLE,
                           total_sdr_sq_err_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
// Get fom tag
std::string fom_tag_name ("fom");
rval = mbi.tag_get_handle(fom_tag_name.c_str(),
                           42,
                           moab::MB_TYPE_DOUBLE,
                           fom_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
std::string total_fom_tag_name ("total_fom");
rval = mbi.tag_get_handle(total_fom_tag_name.c_str(),
                           1,
                           moab::MB_TYPE_DOUBLE,
                           total_fom_tag,
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

// create final meshset
moab::EntityHandle result_meshset;
rval = mbi.create_meshset(moab::MESHSET_SET, result_meshset); MB_CHK_ERR(rval);
MB_CHK_SET_ERR(rval, "Error creating meshset.");


double sdr_contributions = 0.0;
double fom_contributions = 0.0;
double sq_err_sdr = 0.0;
moab::EntityHandle tet_id;
moab::EntityHandle p_src_eh;
double vol;
double adj_flux_h;
double p_src_h;
double sq_err_p_src_h;
double abs_err_p_src_h; 
std::map<int, Tet_info>::iterator mit;
double t_proc = atof(argv[3]);
std::vector<double> sq_err_sdr_vec;//(p_src_map.size(), 0);
std::map<int,double> sq_err_sdr_map;


// for each tet ID, keep running total of the sdr and error
for(mit = p_src_map.begin(); mit!=p_src_map.end(); ++mit){
  tet_id = mit->first;
  vol = p_src_map[tet_id].vol;
  p_src_eh = p_src_map[tet_id].eh;

  std::vector<double> sdr_result(num_e_groups_src, 0);
  std::vector<double> fom_result(num_e_groups_src, 0);
  std::vector<double> sq_err_sdr_result(num_e_groups_src, 0);
  std::vector<double> rel_err_sdr_result(num_e_groups_src, 0);
  std::vector<double> rel_err_p_src_result(num_e_groups_src, 0);
  double sdr_total_result = 0;
  double fom_total_result = 0;
  double sq_err_p_src_total_result = 0;
  double sq_err_sdr_total_result = 0;

  for(int h=0; h <= num_e_groups_src-1; h++){
    adj_flux_h = adj_flux_map[tet_id].result[h];
    p_src_h = p_src_map[tet_id].result[h];
    sq_err_p_src_h = sq_err_p_src_map[tet_id].result[h];
    abs_err_p_src_h = sqrt(sq_err_p_src_h);

    // relative error in photon source per energy group
    if (p_src_h == 0){
      rel_err_p_src_result[h] = 0.0;
    }
    else{
      rel_err_p_src_result[h] = abs_err_p_src_h/p_src_h;
    }
    // accumulate abs error in p src
    sq_err_p_src_total_result += pow(abs_err_p_src_h, 2);

    // sdr contribution per energy group
    sdr_result[h] = adj_flux_h*p_src_h*vol;

    // accumulate sdr contributions over all energy groups
    sdr_total_result += sdr_result[h];

    // sq. abs error in sdr contribution per energy group
    sq_err_sdr_result[h] = pow((adj_flux_h*abs_err_p_src_h*vol), 2);
    // accumulate sq. abs error in sdr over all energy groups
    sq_err_sdr_total_result += sq_err_sdr_result[h];
    if (sdr_result[h] == 0){
      rel_err_sdr_result[h] = 0.0;
      fom_result[h] = 0.0;
    }
    else{
      rel_err_sdr_result[h] = (sqrt(sq_err_sdr_result[h]))/sdr_result[h];
      //fom_result[h] = 1/(sq_err_sdr_result[h]*t_proc);
      fom_result[h] = 1/(pow(rel_err_sdr_result[h] , 2)*t_proc);
    }
  }
  // total abs. err in sdr per tet
  double total_sdr_abs_err = sqrt(sq_err_sdr_total_result);
  // total rel. err in sdr per tet
  double total_sdr_rel_err;
  if(sdr_total_result == 0){ 
    total_sdr_rel_err = 0;
    fom_total_result = 0;
  }
  else{
    total_sdr_rel_err = total_sdr_abs_err/sdr_total_result;
    fom_total_result = 1/( pow(total_sdr_rel_err, 2) * t_proc);
  }

  // accumulate sdr contributions over all energy groups and all  tets
  sdr_contributions += sdr_total_result;
  // accumulate sq. error in sdr over all energy groups and all tets
  sq_err_sdr += sq_err_sdr_total_result;
  sq_err_sdr_vec.push_back(sq_err_sdr_total_result);
  sq_err_sdr_map[tet_id]=sq_err_sdr_total_result;
//  std::cout << "tet id, sq err in sdr " << tet_id << " " << sq_err_sdr_total_result << std::endl;
  // accumulate fom over all tets
  fom_contributions += fom_total_result;

  //set the tag vals 
//  rval = mbi.tag_set_data(p_src_err_tag, &(p_src_eh), 1, &(rel_err_p_src_result[0]));
//  MB_CHK_SET_ERR(rval, "Error setting total err in p src tag val.");
//  double total_psrc_err = sqrt(sq_err_p_src_total_result);
//  rval = mbi.tag_set_data(total_p_src_err_tag, &(p_src_eh), 1, &(total_psrc_err));
//  MB_CHK_SET_ERR(rval, "Error setting total err in p src val.");

//  rval = mbi.tag_set_data(sdr_tag, &(p_src_eh), 1, &(sdr_result[0]));
//  MB_CHK_SET_ERR(rval, "Error setting sdr contribution tag val.");
  rval = mbi.tag_set_data(total_sdr_tag, &(p_src_eh), 1, &(sdr_total_result));
  MB_CHK_SET_ERR(rval, "Error setting total sdr contribution tag val.");

//  rval = mbi.tag_set_data(sdr_err_tag, &(p_src_eh), 1, &(rel_err_sdr_result[0]));
//  MB_CHK_SET_ERR(rval, "Error setting err in sdr contribution tag val.");
  rval = mbi.tag_set_data(total_sdr_err_tag, &(p_src_eh), 1, &(total_sdr_rel_err));
  MB_CHK_SET_ERR(rval, "Error setting err in sdr contribution tag val.");
  rval = mbi.tag_set_data(total_sdr_abs_err_tag, &(p_src_eh), 1, &(total_sdr_abs_err));
  MB_CHK_SET_ERR(rval, "Error setting err in sdr contribution tag val.");
  rval = mbi.tag_set_data(total_sdr_sq_err_tag, &(p_src_eh), 1, &(sq_err_sdr_total_result));
  MB_CHK_SET_ERR(rval, "Error setting err in sdr contribution tag val.");

//  rval = mbi.tag_set_data(fom_tag, &(p_src_eh), 1, &(fom_result[0]));
//  MB_CHK_SET_ERR(rval, "Error setting FOM tag val.");

  rval = mbi.tag_set_data(total_fom_tag, &(p_src_eh), 1, &(fom_total_result));
  MB_CHK_SET_ERR(rval, "Error setting total FOM tag val.");

  //add tet to meshset
  rval = mbi.add_entities(result_meshset, &(p_src_eh), 1);
  MB_CHK_SET_ERR(rval, "Error adding EH to result meshset.");
}

moab::EntityHandle output_list[] = {result_meshset};
rval = mbi.write_mesh("sdr_err_mesh.h5m", output_list, 1);
MB_CHK_SET_ERR(rval, "Error writing mesh.");

double sdr_at_detector = sdr_contributions;
double sdr_abs_err = sqrt(sq_err_sdr);
double sdr_rel_err = sdr_abs_err/sdr_at_detector;
double fom = 1/(pow(sdr_rel_err,2)*t_proc);

//std::vector<double> sorted_sq_err_vec = std::sort(sq_err_sdr_vec.begin(), sq_err_sdr_vec.end());
std::sort(sq_err_sdr_vec.begin(), sq_err_sdr_vec.end());
std::vector<double> cdf; 
std::vector<double>::iterator it;
double cumulative_val; 
std::ofstream cdf_file;
std::ofstream pdf_file;
std::string method = argv[4];
std::string str_tproc = argv[3];
cdf_file.open((method+str_tproc+"_cdf.txt").c_str());
pdf_file.open((method+str_tproc+"_pdf.txt").c_str());
//pdf_file.open((method+std::to_string(t_proc)+"_pdf.txt").c_str());
//for(it = sorted_sq_err_vec.begin(); it !=sorted_sq_err_vec.end(); ++it){
for(it = sq_err_sdr_vec.begin(); it !=sq_err_sdr_vec.end(); ++it){
  cumulative_val += *it;
//  cdf.push_back(cumulative_val);
//  pdf.push_back(*it);
  cdf_file << cumulative_val;
  cdf_file << "\n";
//  std::cout << *it << std::endl;
  pdf_file << *it;
  pdf_file << "\n";
}
cdf_file.close();  
pdf_file.close();  

double sq_err_max = sq_err_sdr_vec.back();
//std::cout << "max err " << sq_err_max << std::endl;
double updated_fom = (pow(sdr_at_detector,2))/(t_proc * (sq_err_sdr-sq_err_max));

std::string ts = argv[2];
//std::cout << ts << " " << sdr_at_detector << " " << sdr_abs_err << " " << sdr_rel_err << " " << fom <<std::endl;
//std::cout << method << " " << t_proc << " " << ts << " " << sdr_at_detector << " " << sdr_abs_err << " " << sdr_rel_err << " " << fom << " " << updated_fom << std::endl;
std::cout << method << " " << t_proc << " " << ts << " " << fom << " " << updated_fom << std::endl;
  
}
