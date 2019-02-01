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
#include "DagMC.hpp"
#include "dagmcmetadata.hpp"

moab::Core mbi;
moab::Tag flux_tag;
moab::Tag error_tag;
moab::Tag id_tag;
moab::Tag move_tag;

moab::DagMC *DAG;
dagmcMetaData* DMD;

struct Tet_info{
       moab::EntityHandle eh;
       std::vector<double> flux;
       std::vector<double> error;
       std::vector<double> sqrd_error;
};

moab::ErrorCode get_mesh_elements(std::string filename,
                                  bool blank,
                                  std::map<int, Tet_info> &tet_map,
                                  moab::EntityHandle &fileset){

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
// Get flux tag
//std::string flux_tag_name ("photon_result");
std::string flux_tag_name ("TALLY_TAG");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//int num_e_groups = num_groups(flux_tag);
int num_e_groups = 24;
std::vector<double> flux(num_e_groups, 0);

//// Get error tag
//std::string error_tag_name ("photon_result_rel_error");
//rval = mbi.tag_get_handle(error_tag_name.c_str(),
//                           moab::MB_TAG_VARLEN,
//                           moab::MB_TYPE_DOUBLE,
//                           error_tag,
//                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
//std::vector<double> error(num_e_groups, 0);
//std::vector<double> sqrd_error(num_e_groups, 0);


// Get all 3D elements in fileset 1
ves.clear();
rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
MB_CHK_SET_ERR(rval, "Error getting tets.");
std::cout << "num ves " << ves.size() << std::endl;
int tet_id;
for (it = ves.begin(); it != ves.end(); ++it){
  if (blank == false){
    // get the flux tag on the ve 
    rval = mbi.tag_get_data(flux_tag, &(*it), 1, &flux[0]);//MB_CHK_ERR(rval);
    MB_CHK_SET_ERR(rval, "Error getting flux tag.");
  }
//  // get the error tag on the ve
//  rval = mbi.tag_get_data(error_tag, &(*it), 1, &error[0] ); MB_CHK_ERR(rval);
//  MB_CHK_SET_ERR(rval, "Error getting error tag.");

  // get the id tag on the ve
  rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id ); MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting id tag.");

  // Add entity handle and flux to tet id map 
  tet_map[tet_id].eh = *it;
  tet_map[tet_id].flux = flux;
  if(blank == false)
    std::cout << "tet id flux" << tet_id << " " << tet_map[tet_id].flux[0] << std::endl;
//  tet_map[tet_id].error = error;
//  tet_map[tet_id].sqrd_error = sqrd_error;
}

 return moab::MB_SUCCESS;
}

int main(int argc, char **argv){

std::string filename1 = argv[1];
std::string filename2 = argv[2];
std::cout << "f name " << filename1 << std::endl;
std::cout << "f name " << filename2 << std::endl;

moab::ErrorCode rval;

moab::EntityHandle fileset;
moab::EntityHandle blankmeshset;

std::map<int, Tet_info> tet_flux_map1;
std::map<int, Tet_info> tet_flux_map2;
std::map<int, Tet_info> diff_flux_map;

int num_e_groups = 24;
bool blank;

  blank = true;
  diff_flux_map.clear();
  rval = get_mesh_elements(filename1, blank, diff_flux_map, blankmeshset);
  MB_CHK_SET_ERR(rval, "Error getting flux mesh file");

  // get flux vals from tet mesh
  blank = false;
  tet_flux_map1.clear();
  rval = get_mesh_elements(filename1, blank, tet_flux_map1, fileset);
  MB_CHK_SET_ERR(rval, "Error getting flux mesh file");
  tet_flux_map2.clear();
  rval = get_mesh_elements(filename2, blank, tet_flux_map2, fileset);
  MB_CHK_SET_ERR(rval, "Error getting flux mesh file");

  // for each tet ID, keep running total of the flux scored in each configuration 
  moab::EntityHandle tet_id;
  std::map<int, Tet_info>::iterator mit;
  for(mit = diff_flux_map.begin(); mit!= diff_flux_map.end(); ++mit){
    tet_id = mit->first;
//    std::cout << "tet1 flux j " << tet_flux_map1[tet_id].flux[0] << std::endl;
//    std::cout << "tet2 flux j " << tet_flux_map2[tet_id].flux[0] << std::endl;
    float avg;
    for(int j=0; j <= num_e_groups-1; j++){
        avg = (tet_flux_map1[tet_id].flux[j]+tet_flux_map2[tet_id].flux[j])/2.0;
        diff_flux_map[tet_id].flux[j] = (tet_flux_map1[tet_id].flux[j]-tet_flux_map2[tet_id].flux[j])/avg;
  //        tot_flux_map[tet_id].sqrd_error[j] += pow(tet_flux_map[tet_id].error[j], 2);
  //        tot_flux_map[tet_id].error[j] = sqrt(tot_flux_map[tet_id].sqrd_error[j]);
    }
//    std::cout << "diff flux j " << diff_flux_map[tet_id].flux[0] << std::endl;
    // set diff flux data tag
    rval = mbi.tag_set_data(flux_tag, &(diff_flux_map[tet_id].eh), 1, &diff_flux_map[tet_id].flux[0]);
    MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
//    rval = mbi.tag_set_data(error_tag, &(tot_flux_map[tet_id].eh), 1, &tot_flux_map[tet_id].error[0]);
//    MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
  }
  // Write out diff flux mesh 
  moab::EntityHandle output_list[] = {blankmeshset};
  rval = mbi.write_mesh("diff_flux.h5m", output_list, 1);
  MB_CHK_SET_ERR(rval, "Error writing out mesh.");

}
