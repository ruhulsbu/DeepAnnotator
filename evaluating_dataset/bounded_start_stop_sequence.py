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

file_list = []
for x in os.listdir('.'):
    if not x.endswith("gz"):
        continue
    file_list.append(x)
    #print x

inframe_write = open("coding_" + sys.argv[1] + "_codon.txt", "w")
outframe_write = open("intragenic_" + sys.argv[1] + "_codon.txt", "w")
generange_write = open("gene_range_" + sys.argv[1] + "_codon.txt", "w")

for file_name in file_list:
    print file_name

    file_read = gzip.open(file_name, "rb")
    
    seq_name = file_read.readline()
    print seq_name
    
    sequence = 'N'*50 + file_read.readline().strip() + 'N'*50
    gene_nogene = '0'*50 + file_read.readline().strip() + '0'*50
    start_codon = '0'*50 + file_read.readline().strip() + '0'*50
    stop_codon = '0'*50 + file_read.readline().strip() + '0'*50

    assert(len(sequence) == len(gene_nogene))
    assert(len(sequence) == len(start_codon))
    assert(len(start_codon) == len(stop_codon))

    """
    for line in file_read:
        print "", len(line)
    #break
    """

    if sys.argv[1] == "start":
        codon_set = ("ATG", "GTG", "TTG", "CTG", "ATT", "ATC")
    else:
        codon_set = ("TAG", "TAA", "TGA") 


    for i in range(50, len(sequence)-52):
        codon = sequence[i-1:i+2]
        #print codon
        
        if codon in codon_set:
            if sys.argv[1] == "start":
                around_gene_boundary = [int(x) for x in start_codon[i-25:i+25]]
            else:
                around_gene_boundary = [int(x) for x in stop_codon[i-25:i+25]]

            if np.sum(around_gene_boundary) == 0:
                continue

            train_seq = sequence[i-50:i+51]
            #print train_seq
            convert = ''.join([seq_convert[train_seq[k]] for k in range(len(train_seq))])

            if gene_nogene[i] == '0':
                 outframe_write.write(convert + '\n')
            else:
                 if sys.argv[1] == "start" and start_codon[i] == '1':
                     generange_write.write(convert + '\n')
                 elif sys.argv[1] == "stop" and stop_codon[i] == '1':
                     generange_write.write(convert + '\n')
                 else:
                     inframe_write.write(convert + '\n')
            """                     
            print convert
            print gene_nogene[i-50:i+51]
            print start_codon[i-50:i+51]
            print stop_codon[i-50:i+51]
            break
            """
    print
    file_read.close()

print "done"
inframe_write.close()
outframe_write.close()
generange_write.close()
