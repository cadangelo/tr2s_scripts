import os
import shutil
import re
from subprocess import call

#src_file = "source_1.h5m"
#src_base = os.path.splitext(src_file)[0]
src_base = "source_1_0"
geom_file = "geom.h5m"
tally_file = "tally.h5m"
trans_file = "tr2s_path.txt"
#main_dir = "/home/tr2s/test_r2s_tally/r2s_test/test_p_setup/nxport"
#lib = "/home/tr2s/basic_mat_lib.h5"
main_dir = "/home/tr2s/test_r2s_tally/r2s_test/test_p_setup/nxport"
pflux_dir = "{0}/pflux".format(main_dir)
lib = "/home/test_r2s_tally/basic_mat_lib.h5"
#map_tag_name = "source_density"
map_tag_name = "photon_result"
blank_tet_mesh = "{0}/box_tetmesh.h5m".format(main_dir) 

os.mkdir(pflux_dir)
for rundir in os.listdir(main_dir):
  if rundir.startswith(src_base):
    os.chdir(rundir)
#    call(["meshtal_to_mesh", "meshtal"])
    #ts = rundir.split('_')[2]
    ts = rundir.split('_')[3]
    #call(["expand_tags.py", "meshtal_tally_4.h5m", "-o",  "pflux_{0}_{1}.vtk".format(src_base, ts)])
    call (["../m2m", "--mesh1", "meshtal_tally_4.h5m", "--mesh2", blank_tet_mesh, "--expand", "false", "--tag", map_tag_name])
    shutil.move("mappeddata.h5m", "{0}/pflux_{1}_{2}.h5m".format(pflux_dir, src_base, ts))
    os.chdir(main_dir)
