import os
import shutil
import re
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call
from matplotlib import pylab

t_proc = str(10000)

tgt_output_file = "tgt_{0}_adj_fom.txt".format(t_proc)
with open(tgt_output_file, 'r') as f:
    data = f.readlines()
tgt_time_step = list()
tgt_fom = list()
tgt_adj_fom = list()
for i, line in enumerate(data):
    meth, tc, ts, fom, adj_fom = line.split()
#    print meth, tc, ts, fom, adj_fom
    tgt_time_step.append(int (ts))
    tgt_fom.append(float (fom))
    tgt_adj_fom.append(float (adj_fom))

c_tgt_output_file = "crop_tgt_{0}_adj_fom.txt".format(t_proc)
with open(c_tgt_output_file, 'r') as f:
    data = f.readlines()
c_tgt_time_step = list()
c_tgt_fom = list()
c_tgt_adj_fom = list()
for i, line in enumerate(data):
    meth, tc, ts, fom, adj_fom = line.split()
#    print meth, tc, ts, fom, adj_fom
    c_tgt_time_step.append(int (ts))
    c_tgt_fom.append(float (fom))
    c_tgt_adj_fom.append(float (adj_fom))

fw_output_file = "fw_{0}_adj_fom.txt".format(t_proc)
with open(fw_output_file, 'r') as f:
    data = f.readlines()
fw_time_step = list()
fw_fom = list()
fw_adj_fom = list()
for i, line in enumerate(data):
    meth, tc, ts, fom, adj_fom = line.split()
#    print meth, tc, ts, fom, adj_fom
    fw_time_step.append(int (ts))
    fw_fom.append(float (fom))
    fw_adj_fom.append(float (adj_fom))

ana_output_file = "ana_{0}_adj_fom.txt".format(t_proc)
with open(ana_output_file, 'r') as f:
    data = f.readlines()
ana_time_step = list()
ana_fom = list()
ana_adj_fom = list()
for i, line in enumerate(data):
    meth, tc, ts, fom, adj_fom = line.split()
#    print meth, tc, ts, fom, adj_fom
    ana_time_step.append(int (ts))
    ana_fom.append(float (fom))
    ana_adj_fom.append(float (adj_fom))


ax = plt.subplot(111)
ax.set_yscale("log")

ax.errorbar(tgt_time_step, tgt_fom, fmt='o', markersize=8, color='royalblue',label="TGT-CADIS")
#ax.errorbar(tgt_time_step, tgt_adj_fom, fmt='*', markersize=10, color='royalblue',label="TGT-CADIS Adjusted")
ax.errorbar(c_tgt_time_step, c_tgt_adj_fom, fmt='*', markersize=10, color='royalblue',label="TGT-CADIS Adjusted")
ax.errorbar(fw_time_step, fw_fom, fmt='o', markersize=8, color='orange', label="FW-CADIS")
ax.errorbar(ana_time_step, ana_fom, fmt='o', markersize=8, color='mediumseagreen', label="Analog")

titname='FOM in Shutdown Dose Rate at Detector over Time';
plt.title(titname,fontsize=14)
plt.xlabel('Time Step',fontsize=18)
plt.ylabel('Figure of Merit [1/min]',fontsize=18)
ax.legend(loc='upper left', shadow=True, ncol=1, numpoints=1)
plt.xlim(-1,17)
plt.xticks(np.arange(0,17,1))

plt.show()
