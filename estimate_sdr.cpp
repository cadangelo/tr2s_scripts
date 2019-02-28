#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <math.h>
#include <stdio.h>
#include <map>

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

struct Tet_info{
       moab::EntityHandle eh;
       std::vector<double> result;
};

int num_groups(moab::Tag tag) {
  moab::ErrorCode rval;
  int tag_size;
  rval = mbi.tag_get_bytes(tag, *(&tag_size));
  if (rval != moab::MB_SUCCESS)
      throw std::runtime_error("Problem getting tag size.");
  return tag_size/sizeof(double);
}

// break string into vector of strings based on delimiter
void tokenize(const std::string& str, 
              std::vector<std::string>& tokens,
              const char* delimiters){
  tokens.clear();

  std::string::size_type next_token_end, next_token_start =
                         str.find_first_not_of( delimiters, 0);

  while ( std::string::npos != next_token_start )
    {
      next_token_end = str.find_first_of( delimiters, next_token_start );
      if ( std::string::npos == next_token_end )
        {
	  tokens.push_back(str.substr(next_token_start));
          next_token_start = std::string::npos;
        }
      else
        {
          tokens.push_back( str.substr( next_token_start, next_token_end -
                                        next_token_start ) );
          next_token_start = str.find_first_not_of( delimiters, next_token_end );
        }
    }
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
rval = mbi.tag_get_handle( GLOBAL_ID_TAG_NAME,
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
rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
MB_CHK_SET_ERR(rval, "Error getting tets.");
std::cout << "Num vol elements: " << ves.size() << std::endl;
int tet_id;
for (it = ves.begin(); it != ves.end(); ++it){
  // get the id tag on the ve
  rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id ); MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting id tag.");

  // get the result tag on the ve 
  rval = mbi.tag_get_data(result_tag, &(*it), 1, &result[0]);
  MB_CHK_SET_ERR(rval, "Error getting result tag.");

  // Add entity handle and result to tet id map 
  tet_map[tet_id].eh = *it;
  tet_map[tet_id].result = result;
}
std::cout << "tet id " << tet_id << std::endl;

return moab::MB_SUCCESS;
}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::EntityHandle fileset;

std::map<int, Tet_info> adj_flux_map;
std::map<int, Tet_info> p_src_map;

// get adjoint photon flux results
std::cout << "Adjoint photon flux file: " << argv[1] << std::endl;
adj_flux_map.clear();
std::string adj_flux_tag_name = "flux";
int num_e_groups_flux;
rval = get_mesh_elements(argv[1], adj_flux_tag_name, adj_flux_map, fileset, num_e_groups_flux);
MB_CHK_SET_ERR(rval, "Error getting adj flux mesh file");

// get forward photon source results
std::cout << "Photon source file: " << argv[2] << std::endl;
p_src_map.clear();
std::string p_src_tag_name = "source_density";
int num_e_groups_src;
rval = get_mesh_elements(argv[2], p_src_tag_name, p_src_map, fileset, num_e_groups_src);
MB_CHK_SET_ERR(rval, "Error getting source mesh file");

//int num_e_groups;
//if(num_e_groups_flux == num_e_groups_src)
//  num_e_groups = num_e_groups_flux;
//else
//  std::cerr << "Number of energy groups from flux file does not match number from source file" << std::endl;

double sdr_at_detector = 0.0;
// for each tet ID, keep running total of the flux scored in each configuration 
moab::EntityHandle tet_id;
std::map<int, Tet_info>::iterator mit;
for(mit = adj_flux_map.begin(); mit!=adj_flux_map.end(); ++mit){
  tet_id = mit->first;
  for(int j=0; j <= num_e_groups_src-1; j++){
      sdr_at_detector += adj_flux_map[tet_id].result[j]*p_src_map[tet_id].result[j];
  }
}

//*** must mult by volume of mesh
std::cout << "SDR at detector = " << sdr_at_detector << std::endl;
  
}
