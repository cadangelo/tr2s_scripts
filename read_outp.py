import os
import shutil
import re
import sys
import matplotlib
import matplotlib.pyplot as plt
from subprocess import call

# outp = sys.argv[1]
# print ('file ', outp)

target = 'cell  14'

#all_dirs = [d[0] for d in os.walk(os.getcwd())]
#source_dirs = list()
#
#for d in all_dirs:
#   last_dir = os.path.split(d)[1]
#   print last_dir
#   if last_dir.startswith('source_1'):
#      source_dirs.append(d)
#
#
## print all_dirs
#print source_dirs

#for source_dir in source_dirs:
results_dir = "results"
os.chdir(results_dir)

time_step = list()
tally = list()
accumulated_dose = list()
error = list()
acc_dose = 0
ts_length = 27.0/32.0

for outp_file in os.listdir(os.getcwd()):
    if outp_file.startswith("outp_"):
        ts = outp_file.split('_')[1] 
        with open(outp_file, 'r') as f:
            data = f.readlines()

        for i, line in enumerate(data):
            if target in line:
                our_data = data[i+24+2]
                energy, tal, err = our_data.split()
                print ts, tal
#                print('time_step:{} energy:{} tally: {} error: {}\n'.format(time_step, energy, tally, error))
                time_step.append(int (ts))
                tally.append(float (tal)*10e-12)
                error.append(float (err))

sorted_tally = [x for _, x in sorted(zip(time_step, tally))]
sorted_time_step = sorted(time_step)
#print sorted_tally
for dr in sorted_tally:
    acc_dose += float(dr)*ts_length
    print dr, acc_dose
    accumulated_dose.append(float (acc_dose))

matplotlib.rc('text',usetex=True)
matplotlib.rc('font',family='serif')

figname='sdr_dose.pdf'
titname='SDR';

# shift the error bars in X
#drp_ebin_diff = 1-np.divide(np.divide(drp_ebin[1]-drp_ebin[0],2),drp_ebin[0])
#drp_ebin_shift = np.multiply(drp_ebin, drp_ebin_diff)

ax = plt.subplot(111)
ax.set_yscale("log")

#drp_b = ax.step(drp_ebin,drp_tal,c='r',label='DRP ON')
#nodrp_b = ax.step(nodrp_ebin,nodrp_tal,c='b',label='DRP ON')

#ax.step(drp_ebin,drp_tal,c='r',label='DRP ON')
#
#ax.errorbar(drp_ebin_shift,drp_tal,yerr=drp_err,fmt=None,ecolor='black')
ax.scatter(time_step, tally)
#ax.scatter(sorted_time_step, accumulated_dose)

plt.xlabel('Time Step',fontsize=18)
plt.ylabel('Shutdown Dose Rate [Sv/s]',fontsize=18)
plt.title(titname,fontsize=14)
plt.xlim(0,32)
#plt.ylim(1e-21,1e-19)

plt.legend()
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
