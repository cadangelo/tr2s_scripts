#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <boost/program_options.hpp>

#include <fstream>
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

namespace po = boost::program_options;

moab::Core mbi;
moab::Tag flux_tag;
moab::Tag id_tag;
moab::Tag move_tag;

moab::DagMC *DAG;
dagmcMetaData* DMD;

struct Tet_info{
       moab::EntityHandle eh;
       std::vector<double> flux;
};

moab::ErrorCode setup(int argc, char* argv[],
                      std::string &tfilename)
//                      std::string &origfilename,
//                      std::vector<std::string> &meshfilenames)
{
  // process command line flags
  po::options_description desc("Allowed options");
  desc.add_options()
      ("help", "Help message")
      ("delta_t", po::value<std::string>(), "Text file with length [s] of each time step")
//      ("mesh", po::value<std::string>(), "Vector of adjoint photon mesh files")
  ;

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);    
  
  if (vm.count("help")) {
      std::cout << desc << std::endl;
      return moab::MB_FAILURE;
  }
  if (vm.count("delta_t")) {
      std::cout << "Delta t file name " << vm["delta_t"].as<std::string>() << std::endl;
      tfilename =  vm["delta_t"].as<std::string>();
  } else {
      std::cout << "Delta t file not given." << std::endl;;
  }
 // if (vm.count("orig")) {
 //     std::cout << "Mesh flux file in original configuration " << vm["orig"].as<std::string>() << std::endl;
 //     origfilename =  vm["orig"].as<std::string>();
 // } else {
 //     std::cout << "Mesh flux file in original configuration not given." << std::endl;;
 // }
 // if (vm.count("mesh")) {
 //     std::cout << "All other mesh flux files " << vm["mesh"].as<std::vector<std::string>>() << std::endl;
 //     meshfilenames =  vm["mesh"].as<std::vector<std::string>>();
 // } else {
 //     std::cout << "Mesh flux files not given." << std::endl;;
 // }

  return moab::MB_SUCCESS;
}

void tokenize( const std::string& str, 
               std::vector<std::string>& tokens,
               const char* delimiters)
{
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

void process_input(std::string tfilename, 
                   std::map< int, double> &ts_delta_map,
                   double &total_time)
{
  std::ifstream transform_input(tfilename);
  std::string line;
  const char* delimiters = " "; 
  const char* time_step_start_token = "s"; 
  int ts;
//  double total_time;

  if (transform_input.is_open())
   {
     while(std::getline(transform_input, line))
       {
         // Skip blank lines in file
         if (line.length() == 0 ) continue;

         // Tokenize the line
         std::vector<std::string> tokens;
         tokenize(line, tokens, delimiters);
         if (tokens.empty()) continue ; 
  
         // Transformation input       
         //if( tokens[0].compare(time_step_start_token ) == 0 && tokens.size() > 1){
             ts = atoi(tokens[0].c_str());
             ts_delta_map[ts] = atof(tokens[1].c_str()); // length of time step
             total_time += ts_delta_map[ts];
         //}
       }
   }
}
moab::ErrorCode sort_filenames(int argc, char **files,
                              std::map<int, std::string> &map_step_to_filename){
  std::string step_num_str;
  int step_num;
  std::string filename;
  std::vector<std::string> tokens;
  const char* delimiter = "_";
  const char* delimiter2 = ".";
  std::vector<std::string>::iterator it;
  int num_files = argc - 3;
  std::cout << "num files from argc "<< num_files << std::endl;
  for (int i = 0 ; i < num_files; ++i){
//  for( it = ofilenames.begin(); it != ofilenames.end(); ++it){
    filename = files[i+3];
    tokenize(filename, tokens, delimiter);
    step_num_str = tokens[2];
    tokenize(step_num_str, tokens, delimiter2);
    step_num = std::atoi(tokens[0].c_str());
    map_step_to_filename[step_num] = filename;
    std::cout << "map s to fname " << step_num << " " << map_step_to_filename[step_num] << std::endl;
  }
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
std::string flux_tag_name ("flux");
rval = mbi.tag_get_handle(flux_tag_name.c_str(),
                           moab::MB_TAG_VARLEN,
                           moab::MB_TYPE_DOUBLE,
                           flux_tag,
                           moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT);
int num_e_groups = num_groups(flux_tag);
std::cout << "num e groups " << num_e_groups << std::endl;
std::vector<double> flux(num_e_groups);


// Get all 3D elements in fileset 1
ves.clear();
rval = mbi.get_entities_by_dimension(fileset, 3, ves);MB_CHK_SET_ERR(rval, "Error getting 3d elements");
MB_CHK_SET_ERR(rval, "Error getting tets.");

int tet_id;
for (it = ves.begin(); it != ves.end(); ++it){
  // get the flux tag on the ve 
  rval = mbi.tag_get_data(flux_tag, &(*it), 1, &flux[0]);//MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting flux tag.");
  // get the id tag on the ve
  rval = mbi.tag_get_data(id_tag, &(*it), 1, &tet_id ); MB_CHK_ERR(rval);
  MB_CHK_SET_ERR(rval, "Error getting id tag.");

  // Add entity handle and flux to tet id map 
  tet_map[tet_id].eh = *it;
  tet_map[tet_id].flux = flux;
}
 std::cout << "flux in last tet id: " << tet_map[tet_id].flux[42] << std::endl;

 return moab::MB_SUCCESS;
}

int main(int argc, char **argv){

moab::ErrorCode rval;

moab::EntityHandle fileset;
moab::EntityHandle blankmeshset;
  
std::string tfilename;
std::string origfilename;
std::vector<std::string> ofilenames;
//rval = setup(argc, argv, tfilename, ofilenames);
rval = setup(argc, argv, tfilename);
MB_CHK_SET_ERR(rval, "Failed to setup meshfile");

double total_time = 0.0;
std::map< int, double> ts_delta_map;
process_input(tfilename, ts_delta_map, total_time);

std::string ofilename;
std::map<int, std::string> step_to_filename;
//sort_filenames(ofilenames, step_to_filename);
sort_filenames(argc, argv, step_to_filename);

std::map<int, Tet_info> tet_flux_map;
std::map<int, Tet_info> tot_flux_map;
//std::map<int, Tet_info> avg_flux_map;

double num_steps = 0.0;
int num_e_groups = 217;
int num_time_steps = ts_delta_map.size();
double delta_t;
std::cout << "tot t, num steps " << total_time << " " << num_time_steps << std::endl;
// for each mesh file on run line
//for (int i = 1; i < argc; ++i){
//for (int ts = 0; ts < num_time_steps; ++ts){
std::map< int, double >::iterator ts_it;
for(ts_it = ts_delta_map.begin(); ts_it != ts_delta_map.end(); ts_it++){
  int ts = ts_it->first;
  delta_t = ts_it->second;
  ofilename = step_to_filename[ts];
  std::cout << "ts, dt, fname " << ts << " " << delta_t << " " << ofilename << std::endl;
//  std::cout << "argc, i " << argc << ", " << i << std::endl; 
//  std::string filename = argv[i];
//  std::cout << filename << std::endl;

  // create tet_id map for blank/final config mesh
  // change this to string comparison for "blankmesh" or "final_config"
  if (ts == 0){
    rval = get_mesh_elements(ofilename, tot_flux_map, blankmeshset);
    MB_CHK_SET_ERR(rval, "Error getting blank mesh file");
    std::cout << "got fname blank map " << ofilename << std::endl;
  }
  else{
    std::cout << "got fname flux map " << ofilename << std::endl;
    // if no time weighting info, will average over number of steps
    //++num_steps;
    // get flux vals from tet mesh
    tet_flux_map.clear();
    rval = get_mesh_elements(ofilename, tet_flux_map, fileset);
    MB_CHK_SET_ERR(rval, "Error getting flux mesh file");
    
    // for each tet ID, keep running total of the flux scored in each configuration 
    std::map<int, Tet_info>::iterator mit;
    moab::EntityHandle tet_id;
    for(mit = tot_flux_map.begin(); mit!=tot_flux_map.end(); ++mit){
      tet_id = mit->first;
      for(int j=0; j <= num_e_groups-1; j++){
        tot_flux_map[tet_id].flux[j] += tet_flux_map[tet_id].flux[j]*delta_t;
        // if this is the last mesh file, find the avg (simplest is just dividing by number of steps)
        if(ts == num_time_steps-1){
            tot_flux_map[tet_id].flux[j] /= total_time;
          }
      }
      if(ts == num_time_steps-1){
            //set the flux tag val to be the avg 
            rval = mbi.tag_set_data(flux_tag, &(tot_flux_map[tet_id].eh), 1, &tot_flux_map[tet_id].flux[0]);
            MB_CHK_SET_ERR(rval, "Error getting flux tag val.");
      }
    }
  }
}

// Write out avg flux mesh 
moab::EntityHandle output_list[] = {blankmeshset};
rval = mbi.write_mesh("avgflux.h5m", output_list, 1);

}
