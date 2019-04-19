import os
import shutil
import re
from subprocess import call

geom_file = "geom_0.h5m"
trans_file = "tr2s_path.txt"
main_dir = "/home/tgt_cadis"
config_file = "config.yml"
lib = "/home/tgt_cadis/basic_mat_lib_2.h5"

call (["./trans", "--geom", geom_file, "--trans", trans_file, "--obbs", "true"])
for filename in os.listdir(main_dir):
  if filename.startswith("geom_"):
    dirName, end = filename.split('.')
    os.mkdir(dirName)
    shutil.move(filename, "{}/geom.h5m".format(dirName))
    shutil.copyfile(config_file, "{}/config.yml".format(dirName))
    os.chdir(dirName)
    call (["/root/opt/dagmc/bin/uwuw_preproc", "-l", lib, "geom.h5m"])
    call(["gtcadis.py", "step1"]) 
    os.chdir(main_dir)

