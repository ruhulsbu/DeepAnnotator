import os, sys, re
import argparse

#parser = argparse.ArgumentParser(description='Parses log file and generates train/val curves')
#parser.add_argument('--log-file', type=str,default="log_tr_va",
#                    help='the path of log file')
#args = parser.parse_args()


TR_RE = re.compile('.*?]\s(1814082, 1) ([\d\.]+)')
VA_RE = re.compile('.*?]\s(604694, 1) ([\d\.]+)')

log = open(sys.argv[2], "r")

log_tr = []
log_va = []
log_tr.append(0.0)
log_va.append(0.0)
for line in log:
    if line.startswith("Training"):
        log_tr.append(float(line.strip().split(" ")[-1]))
    if line.startswith("Evaluating"):
        log_va.append(float(line.strip().split(" ")[-1]))

print_per_epoch = 1
print(len(log_tr), print_per_epoch)
print(log_tr)
print(log_va)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import seaborn as sns
sns.set()

idx = np.arange(len(log_tr))

plt.figure(figsize=(8, 6))
plt.xlabel("Epoch")
plt.ylabel("Accuracy")
plt.plot(idx, log_tr, 'o', linestyle='-', color="r",
         label="Train accuracy")

plt.plot(idx, log_va, 'o', linestyle='-', color="b",
         label="Validation accuracy")
plt.ylim(0.0, 1.01)
plt.legend(loc="best")
#plt.xticks(np.arange(min(idx), max(idx)+1, 1))
#plt.yticks(np.arange(0, 1, 0.2))
#plt.ylim([0,1])
plt.savefig(sys.argv[1]+"_codon_general_train_plot.png")

