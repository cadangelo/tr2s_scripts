import os
import shutil
import math 
import re
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
plt.style.use('seaborn-deep')
#import seaborn as sns
from scipy import stats
from scipy.stats import norm
from subprocess import call
from matplotlib import pylab
#from matplotlib.font_manager import FrontProperties

proc_time = 10000

tgt_time_step = list()
tgt_sdr_t = list()
tgt_error = list()
tgt_rel_error = list()
tgt_cdf = list()
#tgt_output_file = "sdr_t_tgt_{0}.txt".format(proc_time)
#tgt_output_file = "sdr_t_tgt_{0}_box.txt".format(proc_time)
tgt_100_output_file = "tgt_cdf_100.txt"
tgt_1000_output_file = "tgt_cdf_1000.txt"
tgt_10000_output_file = "tgt_cdf_10000.txt"
tgt_100 = list()

with open(tgt_100_output_file, 'r') as f:
    data_100 = f.readlines()
    myint = 100
#    for i, line in enumerate(data_100):
#        tgt_100.append(float(line))
#    	tgt100_100 = [x/myint for x in tgt_100]
with open(tgt_1000_output_file, 'r') as f:
    data_1000 = f.readlines()
with open(tgt_10000_output_file, 'r') as f:
    data_10000 = f.readlines()
    ints = range(0,len(data_10000))



tgt_sig_sq_100 = list()
tgt_sig_sq_1000 = list()
tgt_sig_sq_10000 = list()
fw_sig_sq_100 = list()
fw_sig_sq_1000 = list()
fw_sig_sq_10000 = list()
ana_sig_sq_100 = list()
ana_sig_sq_1000 = list()
ana_sig_sq_10000 = list()

tgt_100_sigma_sq = "tgt_pdf_100.txt"
tgt_1000_sigma_sq = "tgt_pdf_1000.txt"
#tgt_10000_sigma_sq = "tgt_pdf_10000.txt"
tgt_10000_sigma_sq = "tgt_10000_pdf.txt"

#tgt_100_sigma_sq = "tgt_pdf_100_box.txt"
#tgt_1000_sigma_sq = "tgt_pdf_1000_box.txt"
#tgt_10000_sigma_sq = "tgt_pdf_10000_box.txt"

#fw_100_sigma_sq = "fw_pdf_100_box.txt"
#fw_1000_sigma_sq = "fw_pdf_1000_box.txt"
#fw_10000_sigma_sq = "fw_pdf_10000_box.txt"
#fw_100_sigma_sq = "fw_100_pdf.txt"
#fw_1000_sigma_sq = "fw_1000_pdf.txt"
fw_10000_sigma_sq = "fw_10000_pdf.txt"

#ana_100_sigma_sq = "ana_100_pdf.txt"
#ana_1000_sigma_sq = "ana_1000_pdf.txt"
ana_10000_sigma_sq = "ana_10000_pdf.txt"

############## sigma^2 plots
## TGT CADIS
with open(tgt_100_sigma_sq, 'r') as f:
    sig_sq_100 = f.readlines()
    for i, line in enumerate(sig_sq_100):
        tgt_sig_sq_100.append(float(line))
    	tgt_sig_sq_100_100 = [x/100 for x in tgt_sig_sq_100]
with open(tgt_1000_sigma_sq, 'r') as f:
    sig_sq_1000 = f.readlines()
    for i, line in enumerate(sig_sq_1000):
        tgt_sig_sq_1000.append(float(line))
with open(tgt_10000_sigma_sq, 'r') as f:
    sig_sq_10000 = f.readlines()
    for i, line in enumerate(sig_sq_10000):
        tgt_sig_sq_10000.append(float(line))
    	tgt_sig_sq_1000_10 = [x/10 for x in tgt_sig_sq_1000]
## FW CADIS
#with open(fw_100_sigma_sq, 'r') as f:
#    sig_sq_100 = f.readlines()
#    for i, line in enumerate(sig_sq_100):
#        fw_sig_sq_100.append(float(line))
##    	fw_sig_sq_100_100 = [x/100 for x in tgt_sig_sq_100]
#with open(fw_1000_sigma_sq, 'r') as f:
#    sig_sq_1000 = f.readlines()
#    for i, line in enumerate(sig_sq_1000):
#        fw_sig_sq_1000.append(float(line))
with open(fw_10000_sigma_sq, 'r') as f:
    sig_sq_10000 = f.readlines()
    for i, line in enumerate(sig_sq_10000):
        fw_sig_sq_10000.append(float(line))
## Ana CADIS
#with open(ana_100_sigma_sq, 'r') as f:
#    sig_sq_100 = f.readlines()
#    for i, line in enumerate(sig_sq_100):
#        ana_sig_sq_100.append(float(line))
##    	ana_sig_sq_100_100 = [x/100 for x in tgt_sig_sq_100]
#with open(ana_1000_sigma_sq, 'r') as f:
#    sig_sq_1000 = f.readlines()
#    for i, line in enumerate(sig_sq_1000):
#        ana_sig_sq_1000.append(float(line))
with open(ana_10000_sigma_sq, 'r') as f:
    sig_sq_10000 = f.readlines()
    for i, line in enumerate(sig_sq_10000):
        ana_sig_sq_10000.append(float(line))

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',family='serif')
ax = plt.subplot(111)
#ax.set_yscale("log")
#plt.ylim(1e-69,1e-50)
#ax.errorbar(ints, data_100, fmt='x',color='blue',markersize=6, label="TGT-CADIS 100")
#ax.errorbar(ints, tgt100_100, fmt='x',color='black',markersize=6, label="(TGT-CADIS 100)/100")
#ax.errorbar(ints, data_1000, fmt='x',color='green',markersize=6, label="TGT-CADIS 1,000")
#ax.errorbar(ints, data_10000, fmt='x',color='red',markersize=6, label="TGT-CADIS 10,000")

#ax.scatter(ints, tgt_sig_sq_100, fmt='x',color='blue',markersize=6, label="TGT-CADIS 100")
#ax.scatter(ints, tgt_sig_sq_1000, fmt='x',color='green',markersize=6, label="TGT-CADIS 1,000")
#ax.scatter(ints, tgt_sig_sq_10000, fmt='x',color='red',markersize=6, label="TGT-CADIS 10,000")

#mybins = np.logspace(-90,-50,200)
#n,bins,patches = plt.hist(tgt_sig_sq_100, bins=mybins,  log=True, histtype='step',label="100 min")
#plt.hist(tgt_sig_sq_1000, bins=mybins, log=True, histtype='step', label="1,000 min")
#plt.hist(tgt_sig_sq_10000, bins=mybins, log=True, histtype='step', color='red',label="10,000 min")
#plt.hist(tgt_sig_sq_100_100, bins=mybins, log=True, color='orange', histtype='step',label="100 min/100")
#plt.hist([tgt_sig_sq_10000, tgt_sig_sq_1000, tgt_sig_sq_100], log=True,bins=100 label=["10,000 min", "1,000 min", "100 min"])

mybins = np.logspace(-95,-49,100)

# FW
#plt.hist(fw_sig_sq_100, bins=mybins, log=True, color='blue', histtype='step',label="100 min")
#plt.hist(fw_sig_sq_1000, bins=mybins, log=True, color='green', histtype='step',label="FW 1,000 min")
plt.hist(fw_sig_sq_10000, bins=mybins, log=True, color='goldenrod', histtype='step',linewidth=2.1,label="FW-CADIS 10,000 min")

# ANA
#plt.hist(ana_sig_sq_100, bins=mybins, log=True, color='blue', histtype='step',label="100 min")
#plt.hist(ana_sig_sq_1000, bins=mybins, log=True, color='green', histtype='step',label="FW 1,000 min")
plt.hist(ana_sig_sq_10000, bins=mybins, log=True, color='mediumseagreen', histtype='step',linewidth=2.1,label="Analog 10,000 min")

# TGT sigma^2
#plt.hist(tgt_sig_sq_100_100, bins=mybins, log=True, color='silver', histtype='step',linewidth=2.1,label="TGT 100 min/100")
#plt.hist(tgt_sig_sq_1000_10, bins=mybins, log=True, color='skyblue', histtype='step',linewidth=2.1,label="TGT 1,000 min/10")
plt.hist(tgt_sig_sq_10000, bins=mybins, log=True, color='steelblue', histtype='step',linewidth=2.2,label="TGT-CADIS 10,000 min")


#plt.plot(tgt_sig_sq_100, norm.pdf(tgt_sig_sq_100))
#ax.legend(loc='upper left', shadow=True, ncol=1, numpoints=1)
pylab.gca().set_xscale("log")
#pylab.gca().set_yscale("log")
ax.legend(loc='upper right', shadow=True, ncol=1, numpoints=1)
titname='Frequency of SDR Error Values';
plt.title(titname,fontsize=14)
plt.xlabel('$\sigma_2$',fontsize=24)
plt.ylabel('Frequency',fontsize=18)
plt.show()
'''
for i, line in enumerate(data):
    cdf = line
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

fw_time_step = list()
fw_sdr_t = list()
fw_error = list()
fw_rel_error = list()
fw_output_file = "sdr_t_adv_fw_{0}.txt".format(proc_time)
#fw_output_file = "sdr_t_adv_fw_{0}_box.txt".format(proc_time)

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
fw_fom = list()
for re in fw_ts_re:
  fom = 1/(re**2 * proc_time)
  fw_fom.append(fom)

ana_time_step = list()
ana_sdr_t = list()
ana_error = list()
ana_rel_error = list()
ana_output_file = "sdr_t_analog_{0}.txt".format(proc_time)
#ana_output_file = "sdr_t_analog_{0}_box.txt".format(proc_time)

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
ana_fom = list()
for re in ana_ts_re:
  fom = 1/(re**2 * proc_time)
  ana_fom.append(fom)

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
titname='Shutdown Dose Rate at Detector over Time';
#titname='FOM in Shutdown Dose Rate at Detector over Time';

#titname='Accumulated Dose at Detector over Time';

# shift the error bars in X
#drp_ebin_diff = 1-np.divide(np.divide(drp_ebin[1]-drp_ebin[0],2),drp_ebin[0])
#drp_ebin_shift = np.multiply(drp_ebin, drp_ebin_diff)

ax = plt.subplot(111)
ax.set_yscale("log")

#ax.scatter(time_step, tally)
#sdr_err = [e*t for e,t in zip(error,sdr_t)]
#ax.scatter(sorted_time_step, accumulated_dose)

ax.errorbar(tgt_time_step, tgt_sdr_t, yerr=tgt_error, fmt='x',ecolor='blue',markersize=8, label="TGT-CADIS")
ax.errorbar(fw_time_step, fw_sdr_t, yerr=fw_error, fmt='x',color='red',markersize=8,ecolor='orange', label="FW-CADIS")
ax.errorbar(ana_time_step, ana_sdr_t, yerr=ana_error, fmt='x',markersize=8,color='green',ecolor='green', label="Analog")
#ax.errorbar(tgt_time_step, tgt_fom, fmt='+', color='blue',label="TGT-CADIS")
#ax.errorbar(fw_time_step, fw_fom, fmt='+',color='orange', label="FW-CADIS")
#ax.errorbar(ana_time_step, ana_fom, fmt='+',color='green', label="Analog")

plt.xlabel('Time Step',fontsize=18)
#plt.ylabel('Figure of Merit [1/min]',fontsize=18)
plt.ylabel('Shutdown Dose Rate [Sv/s]',fontsize=18)
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
'''
with open(outp) as f:
    for i, line in enumerate(f):
      if target in line:
         target_index = i+1
'''
