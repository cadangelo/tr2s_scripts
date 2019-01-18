import os
import shutil
import re
from subprocess import call

geom_file = "geom_0.h5m"
geom_base = "geom_0"
trans_file = "tr2s_path.txt"
main_dir = "/home/tr2s/gt_cadis"
config_file = "config.yml"
lib = "/home/tr2s/basic_mat_lib.h5"
tetmesh = "box_tetmesh.h5m"
adj_pflux_dir = "{}/adj_pflux".format(main_dir)
adj_v = "adj_p.h5m"
t_matrix = "step2_T.npy"

os.mkdir(adj_pflux_dir)
for rundir in os.listdir(main_dir):
  if rundir.startswith(geom_base):
    shutil.copyfile(tetmesh, "{}/tetmesh.h5m".format(rundir))
    shutil.copyfile(tetmesh, "{0}/{1}".format(rundir, t_matrix))
    os.chdir(rundir)
    ts = rundir.split('_')[2]
    call(["python", "../atflux_to_mesh.py"])
    call(["../m2m", "--expand", "false", "--mesh1", adj_v, "--mesh2", "tetmesh.h5m", "--tag", "flux"]) 
    shutil.move("mappeddata.h5m", "{0}/adj_pflux_{1}_{2}.h5m".format(adj_pflux_dir, geom_base, ts))
    os.chdir(main_dir)

