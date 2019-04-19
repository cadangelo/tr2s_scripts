import os
import shutil
import re
from subprocess import call

geom_file = "geom_0.h5m"
#geom_base = "geom_0"
geom_base = "nw_geom_0"
trans_file = "tr2s_path.txt"
main_dir = "/home/tgt_cadis"
config_file = "config.yml"
lib = "/home/tgt_cadis/basic_mat_lib.h5"
adj_pflux_dir = "{}/adj_pflux".format(main_dir)
adj_v = "adj_p.h5m"
t_matrix = "step2_T.npy"
mesh_file = "tr2s_unstr_scaled_better.h5m"
mesh_base = mesh_file.split('.')[0]

#os.mkdir(adj_pflux_dir)
#shutil.copyfile(mesh_file, "{0}_0.h5m".format(mesh_base))
#call (["./trans", "--geom", mesh_file, "--trans", trans_file, "--obbs", "false"])
for rundir in os.listdir(main_dir):
  if rundir.startswith(geom_base):
    #ts = rundir.split('_')[2]
    ts = rundir.split('_')[3]
    shutil.move("{0}_{1}.h5m".format(mesh_base, ts), "{0}/tetmesh.h5m".format(rundir))
    shutil.copyfile(t_matrix, "{0}/{1}".format(rundir, t_matrix))
    os.chdir(rundir)
    print ("current run dir ", rundir)
    call(["python", "../atflux_to_mesh.py"])
    call(["../m2m", "--expand", "false", "--mesh1", adj_v, "--mesh2", "tetmesh.h5m", "--tag", "flux"]) 
    shutil.copyfile("mappeddata.h5m", "{0}/adj_pflux_{1}.h5m".format(adj_pflux_dir, ts))
    os.chdir(main_dir)

