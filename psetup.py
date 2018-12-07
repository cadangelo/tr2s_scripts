import os
import shutil
import re
from subprocess import call

src_file = "source_1.h5m"
src_base, ext = os.path.splitext(src_file)
geom_file = "geom.h5m"
tally_file = "tally.h5m"
trans_file = "tr2s_path.txt"
p_inp_file = "photon_transport.inp"
ebounds_file = "e_bounds"
main_dir = "/home/tr2s/fw_p/tr_srcs_2/test_p_setup"

call (["./trans", src_file, trans_file])
for filename in os.listdir(main_dir):
  if filename.startswith("source_"):
    dirName, end = filename.split('.')
    os.mkdir(dirName)
    shutil.move(filename, "{}/source.h5m".format(dirName))
    shutil.copyfile(p_inp_file, "{}/photon_dag.inp".format(dirName))
    shutil.copyfile(ebounds_file, "{}/e_bounds".format(dirName))

#    os.chdir(dirName)
#    print (os.getcwd())
    

call (["./trans", geom_file, trans_file])
call (["./trans", tally_file, trans_file])
for filename in os.listdir(main_dir):
  if filename.startswith("geom_", "tally_"):
    start, ts = (os.path.splitext(filename))[0].split('_')
    shutil.move(filename, "{0}_{1}/{2}.h5m".format(src_base, ts, start))
