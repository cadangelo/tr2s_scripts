import os
import shutil
import re
from subprocess import call

main_dir = "/home/tgt_cadis/split_chamber/tgt_vr"
p_src_file = "{0}/fw_n/source_1.h5m".format(main_dir)
p_src_base = os.path.splitext(os.path.basename(p_src_file))[0]
sq_err_file = "{0}/fw_n/sq_err_p_src_m.h5m".format(main_dir)
sq_err_base = os.path.splitext(os.path.basename(sq_err_file))[0]
trans_file = "{0}/tr2s_path.txt".format(main_dir)
geom_file = "geom_0.h5m"
geom_base = "geom_0"
adj_p_f_mesh = "adj_p.h5m"
sdr_output = "{0}/sdr_t_output.txt".format(main_dir)
global_mesh = "{0}/fw_n/blank_mesh.h5m".format(main_dir)

#call (["./trans", "--geom", p_src_file, "--trans", trans_file, "--obbs", "false"])
#call (["./trans", "--geom", sq_err_file, "--trans", trans_file, "--obbs", "false"])
#shutil.copyfile(p_src_file, "{0}_0.h5m".format(p_src_base))
shutil.copyfile(sq_err_file, "{0}_0.h5m".format(p_src_base))
with open(sdr_output, "a") as output:
  for rundir in os.listdir(main_dir):
    if rundir.startswith(geom_base):
      ts = rundir.split('_')[2]
      shutil.copyfile("{0}_{1}.h5m".format(p_src_base, ts), "{0}/psrc_mesh.h5m".format(rundir))
      shutil.copyfile("{0}_{1}.h5m".format(sq_err_base, ts), "{0}/err_mesh.h5m".format(rundir))
      os.chdir(rundir)
      call(["../m2m", "--expand", "false", "--mesh1", "adj_p.h5m", "--mesh2", "../global_vox_mesh.h5m", "--tag", "flux"]) 
      #call(["../m2m", "--expand", "false", "--mesh1", "adj_p.h5m", "--mesh2", global_mesh, "--tag", "flux"]) 
      shutil.move("mappeddata.h5m", "adj_p_vox.h5m")
      call(["../m2m", "--expand", "false", "--mesh1", "psrc_mesh.h5m", "--mesh2", "adj_p_vox.h5m", "--tag", "source_density"]) 
      shutil.move("mappeddata.h5m", "adj_p_fw_src_vox.h5m")
      call(["../m2m", "--expand", "false", "--mesh1", "err_mesh.h5m", "--mesh2", "adj_p_fw_src_vox.h5m", "--tag", "sq_err_q_src"]) 
      shutil.move("mappeddata.h5m", "adj_p_fw_src_err_vox.h5m")
      #outfile.write("Time step: %s " % ts)
      call(["../sdr", "adj_p_fw_src_err_vox.h5m", ts], stdout=output)
      os.chdir(main_dir)
output.close()
