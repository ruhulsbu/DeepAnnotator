# -*- coding: utf-8 -*-
"""
Created on Mon Feb  5 14:29:25 2018

@author: moamin
"""
import sys, os
import numpy as np
import gzip

seq_convert = {'N':'0', 'A':'1', 'C':'2', 'G':'3', 'T':'4'}
print seq_convert
maxlen = 101

file_list = []
for x in os.listdir('.'):
    if not (x.startswith("batch") and x.endswith("gz")):
        continue
    file_list.append(x)
    #print x

genome_write = gzip.open("whole_genome_sequences.gz", "w")

for file_name in file_list:
    print file_name

    file_read = gzip.open(file_name, "rb")
    
    seq_name = file_read.readline()
    sequence = 'N'*maxlen + file_read.readline().strip() + 'N'*maxlen
    geneinfo = '0'*maxlen + file_read.readline().strip() + '0'*maxlen

    print seq_name, len(sequence) 

    assert(len(sequence) == len(geneinfo))
    genome_write.write(">"+file_name+"\n")

 
    for i in range(0, len(sequence)-maxlen+1):
        sub_seq = sequence[i:i+maxlen]
        convert = ''.join([seq_convert[sub_seq[k]] for k in range(len(sub_seq))])
        
        sub_seq = geneinfo[i:i+maxlen]
        zeros = sub_seq.count('0')
        
        genome_write.write(convert+'\n')

    print i
    file_read.close()

    #break

print "done"
genome_write.close()

