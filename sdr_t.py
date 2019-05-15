import os
import shutil
import re
from subprocess import call

main_dir = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/tgt_vr"
t_proc = str(10000)
# tgt == 1; fw == 2; analog == 3
method = 'ana' 
if (method == 'tgt'):
	## 2 step
	#p_src_err_file = "{0}/test_2_steps/fw_n/{1}/sq_err_p_src.h5m".format(main_dir, t_proc)
	#p_src_err_file = "{0}/test_2_steps/fw_n/{1}/sq_err_p_src_box.h5m".format(main_dir, t_proc)
	### TGT
	#p_src_err_file = "{0}/fw_n/{1}/sq_err_p_src_box.h5m".format(main_dir, t_proc)
	p_src_err_file = "{0}/fw_n/{1}/sq_err_p_src.h5m".format(main_dir, t_proc)
if (method == 'fw'):
	### FW
	#p_src_err_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/fw_vr/adv/{1}/sq_err_p_src_box.h5m".format(main_dir, t_proc)
	p_src_err_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/fw_vr/adv/{1}/sq_err_p_src.h5m".format(main_dir, t_proc)
if (method == 'ana'):
	### ANALOG
	#p_src_err_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/analog/{1}/sq_err_p_src_box.h5m".format(main_dir, t_proc)
	p_src_err_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/analog/{1}/sq_err_p_src.h5m".format(main_dir, t_proc)

p_src_err_base = os.path.splitext(os.path.basename(p_src_err_file))[0]
print("file name ", p_src_err_file)
#p_src_file = "{0}/fw_n/10000/source_1.h5m".format(main_dir)
#p_src_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/fw_vr/10000/source_1.h5m".format(main_dir)
#p_src_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/analog/10000/source_1.h5m".format(main_dir)
#p_src_base = os.path.splitext(os.path.basename(p_src_file))[0]
#sq_err_file = "{0}/fw_n/10000/sq_err_p_src.h5m".format(main_dir)
#sq_err_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/fw_vr/10000/sq_err_p_src.h5m".format(main_dir)
#sq_err_file = "/home/tgt_cadis/split_chamber/far_det_thicc_wall/analog/10000/sq_err_p_src.h5m".format(main_dir)
#sq_err_base = os.path.splitext(os.path.basename(sq_err_file))[0]


#*** change trans file!!
trans_file = "{0}/tr2s_path.txt".format(main_dir)
#trans_file = "two_step.txt"
geom_file = "geom_0.h5m"
geom_base = "geom_"
adj_p_f_mesh = "adj_p.h5m"
#*** change output file
#sdr_output = "{0}/sdr_t_output.txt".format(main_dir)
sdr_output = "{0}/{1}_{2}_adj_fom.txt".format(main_dir,method,t_proc)
global_mesh = "{0}/fw_n/blank_mesh.h5m".format(main_dir)

vr_folder = "{0}/del2_{1}_srcs_err_{2}".format(main_dir,method, t_proc)
vtk_folder = "{0}/del_tgt_{1}_vtk".format(main_dir, t_proc)
#vr_folder = "{0}/terr_analog_srcs_err_{1}".format(main_dir, t_proc)
#vtk_folder = "{0}/terr_analog_{1}_vtk".format(main_dir, t_proc)
#vtk_folder = "{0}/te3_adv_fw_{1}_vtk".format(main_dir, t_proc)
#vr_folder = "{0}/te3_adv_fw_srcs_err_{1}".format(main_dir, t_proc)


#call(["/home/tgt_cadis/m2m", "--expand", "false", "--mesh1", sq_err_file, "--mesh2", p_src_file, "--tag", "sq_err_p_src"])
#shutil.move("mappeddata.h5m", "{0}.h5m".format(sq_err_base)) 
#call (["/home/tgt_cadis/trans", "--geom", p_src_file, "--trans", trans_file, "--obbs", "false"])
#call (["/home/tgt_cadis/trans", "--geom", sq_err_file, "--trans", trans_file, "--obbs", "false"])

call (["/home/tgt_cadis/trans", "--geom", p_src_err_file, "--trans", trans_file, "--obbs", "false"])
shutil.copyfile(p_src_err_file, "{0}_0.h5m".format(p_src_err_base))
os.mkdir(vr_folder)
for filename in os.listdir(main_dir):
  if filename.startswith(p_src_err_base):
    shutil.move(filename, "{0}/{1}".format(vr_folder, filename))
#os.mkdir(vtk_folder)

with open(sdr_output, "a") as output:
  for rundir in os.listdir(main_dir):
    if rundir.startswith(geom_base):
    #if rundir.startswith("geom_10"):
      ts = rundir.split('_')[1]
      print ("ts ", ts)
      os.chdir(rundir)
      ### change!
      call(["/home/tgt_cadis/m2m", "--expand", "false", "--mesh1", "{0}/sq_err_p_src_{1}.h5m".format(vr_folder, ts), "--mesh2", "adj_p_flux.h5m", "--tag", "est_p_src"]) 
      #call(["/home/tgt_cadis/m2m", "--expand", "false", "--mesh1", "{0}/sq_err_p_src_1.h5m".format(vr_folder), "--mesh2", "adj_p_flux.h5m", "--tag", "est_p_src"]) 
      shutil.move("mappeddata.h5m", "adj_p_flux_src.h5m")
      ### change!
      call(["/home/tgt_cadis/m2m", "--expand", "false", "--mesh1", "{0}/sq_err_p_src_{1}.h5m".format(vr_folder, ts), "--mesh2", "adj_p_flux_src.h5m", "--tag", "sq_err_p_src"]) 
      #call(["/home/tgt_cadis/m2m", "--expand", "false", "--mesh1", "{0}/sq_err_p_src_1.h5m".format(vr_folder), "--mesh2", "adj_p_flux_src.h5m", "--tag", "sq_err_p_src"]) 
      shutil.move("mappeddata.h5m", "adj_p_flux_src_err.h5m")
##      call(["/home/tgt_cadis/sdr", "adj_p_flux_src_err.h5m", ts, t_proc], stdout=output)
#      shutil.move("sdr_err_mesh.h5m", "tbox_sdr_mesh.h5m")
#      call(["/root/opt/moab/bin/mbconvert", "sdr_err_mesh.h5m", "{0}/tot_sdr_{1}.vtk".format(vtk_folder,ts) ])
#      shutil.move("sdr_err_mesh.h5m", "total_sdr_mesh.h5m")
      call(["/home/tgt_cadis/fom", "adj_p_flux_src_err.h5m", ts, t_proc,method], stdout=output)
#      call(["expand_tags.py", "-t", "sdr", "sdr_err_mesh.h5m", "-o", "{0}/sdr_{1}.vtk".format(vtk_folder, ts) ])
#      call(["expand_tags.py", "-t", "sdr_err", "sdr_err_mesh.h5m", "-o", "{0}/sdr_err_{1}.vtk".format(vtk_folder,ts) ])
#      call(["expand_tags.py", "-t", "p_src_err", "sdr_err_mesh.h5m", "-o", "{0}/p_src_err_{1}.vtk".format(vtk_folder,ts) ])
#      call(["expand_tags.py", "-t", "est_p_src", "sdr_err_mesh.h5m", "-o", "{0}/est_p_src_{1}.vtk".format(vtk_folder,ts) ])
#      call(["expand_tags.py", "-t", "fom", "sdr_err_mesh.h5m", "-o", "{0}/fom_{1}.vtk".format(vtk_folder,ts) ])
#      call(["/root/opt/moab/bin/mbconvert", "sdr_err_mesh.h5m", "{0}/sq_abs_{1}.vtk".format(vtk_folder,ts) ])
      os.chdir(main_dir)
output.close()
