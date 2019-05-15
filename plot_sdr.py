import os
import shutil
import re
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from matplotlib import pylab
#from matplotlib.font_manager import FrontProperties

proc_time = 10000

tgt_time_step = list()
tgt_sdr_t = list()
tgt_error = list()
tgt_rel_error = list()
#tgt_output_file = "sdr_t_tgt_{0}.txt".format(proc_time)
tgt_output_file = "sdr_t_tgt_{0}_box.txt".format(proc_time)

with open(tgt_output_file, 'r') as f:
    data = f.readlines()

for i, line in enumerate(data):
    ts, sdr, err, rel_err = line.split()
    print ts, sdr 
    #print('time_step:{} energy:{} tally: {} error: {}\n'.format(time_step, energy, tally, error))
    tgt_time_step.append(int (ts))
    # orig tally units are Sv/s-- convert to microSv/hr
    tgt_sdr_t.append(float (sdr))
    tgt_error.append(float (err))
    tgt_rel_error.append(float (rel_err))

tgt_ts_re = [x for _, x in zip(tgt_time_step, tgt_rel_error)]
tgt_fom = list()
for re in tgt_ts_re:
  fom = 1/(re**2 * proc_time)
  tgt_fom.append(fom)
  print ("tgt fom", fom)

fw_time_step = list()
fw_sdr_t = list()
fw_error = list()
fw_rel_error = list()
#fw_output_file = "sdr_t_adv_fw_{0}.txt".format(proc_time)
fw_output_file = "sdr_t_adv_fw_{0}_box.txt".format(proc_time)

with open(fw_output_file, 'r') as f:
    data = f.readlines()

for i, line in enumerate(data):
    ts, sdr, err, rel_err = line.split()
    print ts, sdr 
    #print('time_step:{} energy:{} tally: {} error: {}\n'.format(time_step, energy, tally, error))
    fw_time_step.append(int (ts))
    # orig tally units are Sv/s-- convert to microSv/hr
    fw_sdr_t.append(float (sdr))
    fw_error.append(float (err))
    fw_rel_error.append(float (rel_err))

fw_ts_re = [x for _, x in zip(fw_time_step, fw_rel_error)]
print("fw ts", fw_time_step)
fw_fom = list()
for re in fw_ts_re:
  fom = 1/(re**2 * proc_time)
  fw_fom.append(fom)
  print ("fw fom", fom)

ana_time_step = list()
ana_sdr_t = list()
ana_error = list()
ana_rel_error = list()
#ana_output_file = "sdr_t_analog_{0}.txt".format(proc_time)
ana_output_file = "sdr_t_analog_{0}_box.txt".format(proc_time)

with open(ana_output_file, 'r') as f:
    data = f.readlines()

for i, line in enumerate(data):
    ts, sdr, err, rel_err = line.split()
    print ts, sdr 
    #print('time_step:{} energy:{} tally: {} error: {}\n'.format(time_step, energy, tally, error))
    ana_time_step.append(int (ts))
    # orig tally units are Sv/s-- convert to microSv/hr
    ana_sdr_t.append(float (sdr))
    ana_error.append(float (err))
    ana_rel_error.append(float (rel_err))
   
ana_ts_re = [x for _, x in zip(ana_time_step, ana_rel_error)]
print("ana ts", ana_time_step)
ana_fom = list()
for re in ana_ts_re:
  fom = 1/(re**2 * proc_time)
  ana_fom.append(fom)
  print ("ana fom", fom)

#sorted_tally = [x for _, x in sorted(zip(time_step, sdr_t))]
#sorted_time_step = sorted(time_step)
#ts_length = 27.0/17.0
#accumulated_dose = list()
#acc_dose = 0
#print ("sorted tal", sorted_tally)
#for dr in sorted_tally:
#    acc_dose += float(dr)*ts_length
#    print ("dr ", dr)
#    accumulated_dose.append(float (acc_dose))

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',family='serif')

#figname='sdr_dose.pdf'
#figname='sdr_tgt_fw_analog.pdf'
figname='sdr_fom_tgt_fw_analog_{0}.pdf'.format(proc_time)
#titname='Shutdown Dose Rate at Detector over Time';
titname='FOM in Shutdown Dose Rate at Detector over Time';

#titname='Accumulated Dose at Detector over Time';

# shift the error bars in X
#drp_ebin_diff = 1-np.divide(np.divide(drp_ebin[1]-drp_ebin[0],2),drp_ebin[0])
#drp_ebin_shift = np.multiply(drp_ebin, drp_ebin_diff)

ax = plt.subplot(111)
ax.set_yscale("log")

#ax.scatter(time_step, tally)
#sdr_err = [e*t for e,t in zip(error,sdr_t)]
#ax.scatter(sorted_time_step, accumulated_dose)

#ax.errorbar(tgt_time_step, tgt_sdr_t, yerr=tgt_error, fmt='x',ecolor='blue',markersize=8, label="TGT-CADIS")
#ax.errorbar(fw_time_step, fw_sdr_t, yerr=fw_error, fmt='x',color='red',markersize=8,ecolor='orange', label="FW-CADIS")
#ax.errorbar(ana_time_step, ana_sdr_t, yerr=ana_error, fmt='x',markersize=8,color='green',ecolor='green', label="Analog")

ax.errorbar(tgt_time_step, tgt_fom, fmt='+', color='blue',label="TGT-CADIS")
ax.errorbar(fw_time_step, fw_fom, fmt='+',color='orange', label="FW-CADIS")
ax.errorbar(ana_time_step, ana_fom, fmt='+',color='green', label="Analog")

plt.xlabel('Time Step',fontsize=18)
plt.ylabel('Figure of Merit [1/min]',fontsize=18)
#plt.ylabel('Shutdown Dose Rate [Sv/s]',fontsize=18)
#plt.ylabel('Dose [microSv]',fontsize=18)
plt.title(titname,fontsize=14)
plt.xlim(-1,17)
plt.xticks(np.arange(0,17,1))

#plt.legend()
#ax.legend(loc='upper center', bbox_to_anchor=(0.3, 1.0),ncol=3)

ax.legend(loc='upper left', shadow=True, ncol=1, numpoints=1)
#ax.legend(loc='best', bbox_to_anchor=(0.3,1.0), shadow=True, ncol=1, numpoints=1)
plt.show()

#pylab.legend(loc=9, bbox_to_anchor=(0.5,-0.1), ncol=3)
#plt.savefig(figname)
plt.show()
#plt.clf()

# print our_data
'''
with open(outp) as f:
    for i, line in enumerate(f):
      if target in line:
         target_index = i+1
'''
