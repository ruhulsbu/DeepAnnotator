from os import walk
import sys, gzip
import csv, re, random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from mpi4py import MPI
#comm = MPI.COMM_WORLD
#name=MPI.Get_processor_name()
#print("hello world")
#print(("name :",name,"my rank is: ",comm.rank,"out of total: ",comm.size))

genmark_organisms = ['Bacillus subtilis', 'Escherichia coli', 'Haemophilus influenzae', \
			'Helicobacter pylori', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae', \
			'Methanobacterium thermoauthotrophicum', 'Synechocystis', \
			'Archeoglobus fulgidus', 'Methanocaldococcus jannaschii']


file_read = open("./dataset/assembly_summary.txt", "r")
complete = {}
unique_species = {}
for line in file_read:
	columns = line.split("\t")
	if len(columns) < 21:
		print line
		continue

	species_prefix = columns[7].split(" ")[0]
	genmark_flag = False
	for k in range(0, len(genmark_organisms)):
		if columns[7].startswith(genmark_organisms[k]):
			genmark_flag = True
			#species_prefix = genmark_organisms[k]
			break
	if genmark_flag == True:
		continue

	if unique_species.has_key(species_prefix):
		#print columns
		continue
	else:
		if columns[11] == "Complete Genome":
			#print columns
			print columns[0], columns[7], columns[11]
			complete[columns[0]] = True
			unique_species[species_prefix] = True
print "Length of Complete = ", len(complete)
#exit()

fna_dir_name = "./dataset/GbBac_FNA"#"debug_files/"#"dataset/GbBac_FNA"#'test_fna'
gff_dir_name = "./dataset/GbBac_GFF"#"debug_files/"#"dataset/GbBac_GFF"#'test_gff'
output_dir_name = "./genmark_train/"#"./batch_files"#"./test_files"

gff_file_list = []
for (dirpath, dirnames, filenames) in walk(gff_dir_name):
	for i in range(0, len(filenames)):
		if complete.has_key("GCA_" + filenames[i].split("_")[1]) == True and "gff" in filenames[i]:
    			gff_file_list.append(filenames[i])
    	break

#gff_file_list = ["GCA_000764535.1_ASM76453v1_genomic.gff.gz"]
#print "GFF File Count: ", len(gff_file_list)
#exit()

# read in the synonymous codon map, no start codons, no stop codons (19 Amino Acids Total)
aa_list = []
codon_to_aa_dict = dict()
csvfile = open('AA_CODONS_MAP_no_stop')
reader = csv.reader(csvfile, delimiter='\t')
for row in reader:
    mylist = row[3].split(',')
    aa_list.append(row[0])
    for codon in mylist:
        codon_to_aa_dict[codon] = row[0]

begind = 0
endind = len(gff_file_list)

comm = MPI.COMM_WORLD
bin_length = 1 + len(gff_file_list) / comm.size
begind = comm.rank * bin_length
endind = (comm.rank + 1) * bin_length
#endind = 100#debug
print "Processing: ", (comm.rank, begind, min(endind, len(gff_file_list)), comm.size)
#exit()


debug = 9
count = 0
wrong = 0

file_write = None
batch_size = 1


################################################################################

def annotate_sequence_following_gff(chrom_letters, csvfile):
	current_chrom = ""
	# dictionary to hold strings of 1's and 0's corresponding to labels for this chromosome/region
	label_dict = dict()
	region_dict = {}

	# process annotation, just for this first region
	#current_chrom = ""
	reader = csv.reader(csvfile, delimiter='\t')
	for row in reader:
    		if (row[0][0] == "#"):
        		continue # a comment line
    
    		if (row[2] == "region"):
        		if (current_chrom != ""): # ****we are done with the whole region - output everything
            			# here, output to file all the labels for this region
            			# TODO:  output just the sequence, or perhaps do the three-mer or hepta-mer conversion here???
            			# TODO: caveat - will not print out the "last" region with this loop; repeat under the loop?
				"""
            			file_write.write(">" + current_chrom + "\n")
            			file_write.write(chrom_letters + "\n")
            			# here is is important that we keep the order of printing out the label_dict consistent
            			# this will be the order in which we will interpret the labels in the verification pipeline
            			file_write.write(label_dict["gene"] + "\n")
            			file_write.write(label_dict["start_codon"] + "\n")
            			file_write.write(label_dict["stop_codon"] + "\n")
            			file_write.write(label_dict["upstream"] + "\n")
            			file_write.write(label_dict["downstream"] + "\n")
            			file_write.write(label_dict["intragenic"] + "\n")
            			file_write.write(label_dict["exon"] + "\n")
            			file_write.write(label_dict["intron"] + "\n")
            			for aa in aa_list: # in sorted alphabetical order
                			file_write.write(label_dict[aa] + "\n")
            			"""
				#label_dict["intragenic"] = label_dict["intragenic"][:prev_gene_stop] + \
				#		"1" * (len(chrom_letters[current_chrom])-prev_gene_stop)

				assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))
				region_dict[current_chrom] = label_dict
				print "Finish Annotating Region: ", (current_chrom, len(label_dict["gene"]), len(chrom_letters[current_chrom]), len(region_dict[current_chrom]["sequence"]))
			
	        	# *** start the new region processing
        		label_dict = {}  # remove all entries from the labels dictionary
        		chrom_length = int(row[4]) - int(row[3]) + 1
        		current_chrom = row[0]
			# make sure the gff region and chrom_letters are the same length
			chrom_length = min(len(chrom_letters[current_chrom]), chrom_length, 102400000)
			if len(chrom_letters[current_chrom]) > chrom_length:
				chrom_letters[current_chrom] = chrom_letters[current_chrom][:chrom_length]
			assert(len(chrom_letters[current_chrom]) == chrom_length)
			label_dict["sequence"] = chrom_letters[current_chrom]			
	
        		# label for gene or not gene
        		label_dict["gene"] = "0" * chrom_length
        		# labels for start codon and stop codon
       		 	label_dict["start_codon"] = "0" * chrom_length
       			label_dict["stop_codon"] = "0" * chrom_length
			"""
        		# labels for upstream, downstream, and intragenic regions
       			label_dict["upstream"] = "0" * chrom_length
        		label_dict["downstream"] = "0" * chrom_length
        		label_dict["intragenic"] = "0" * chrom_length
        		# labels for coding and non-coding parts of the gene
        		# for now, call them exon and intron
        		label_dict["exon"] = "0" * chrom_length
        		label_dict["intron"] = "0" * chrom_length
        		label_dict["frame1"] = "0" * chrom_length
        		# also, make label strings for each Amino Acid
        		# make a label string for this amino_acid
        		for aa in aa_list:
            			label_dict[aa] = "0" * chrom_length
			"""
        		# reset some other variables
        		prev_gene_stop = 0 # used to set the intragenic regions
        		consider_coding_regions = 1 # set to 0 when we read coding regions of the gene that we skipped
       			prev_coding_start_end = (0, 0) # used to set the intron regions and sanity checks
			assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))

		if (row[2] == "gene"):
		
			# Important: for now, only consider protein-coding genes 
			# - tRNA genes etc do not start with ATG, or end with stop
			# TODO:  do sanity check here - should be keeping the majority of genes
			if (("gene_biotype=protein_coding" not in row[8]) or (row[6] == "-")):
				#print("skipping: ", row)
				consider_coding_regions = 0 # do not consider exons or CDSs for this gene!!!
				continue
		
			consider_coding_regions = 1
		
			#print(row)
			gene_start = int(row[3])-1 # subtract 1 because annotation uses 1-based coordinates
			gene_end = int(row[4])
			gene_length = gene_end - gene_start

			#if (gene_start > 30000):#Condition for debugging
			#	break

			if (gene_end > chrom_length):
				continue
		
			# assign the gene label to 1's for the stretch of this gene
			#print(label_dict["gene"][gene_start-1:gene_end+1])
			label_dict["gene"] = label_dict["gene"][:gene_start] + "1" * gene_length \
				    + label_dict["gene"][gene_start+gene_length:]
			assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))
			#print(label_dict["gene"][gene_start-1:gene_end+1])
		
			# for now, assume that the first 3 nt of a gene are ATG (or GTG or TTG)
			# and the last 3 are STOP, still label but output error otherwise to debug
			# check that the start codon is the beginning of the gene and set the label
			if (chrom_letters[current_chrom][gene_start:gene_start+3] in ("ATG", "GTG", "TTG", "CTG", "ATT", "ATC")):
				label_dict["start_codon"] = label_dict["start_codon"][:gene_start]  + "111" + \
						label_dict["start_codon"][gene_start+3:]
				#print(label_dict["start_codon"][gene_start-1:gene_end+1])
				#print(row)
			else:
				# NOTE:  I see a few of these from time to time, I don't think this is a problem necessarily
				# as long as the stop codon is correct - that one is more important
				#print("Possible problem: No start codon at the start of the gene, not labeling\n", row) #TODO:  add file name to error message
				pass
			# check that the stop codon is the end of the gene and set the label
			if (chrom_letters[current_chrom][gene_end-3:gene_end] in ("TAG", "TAA", "TGA")) :
				label_dict["stop_codon"] = label_dict["stop_codon"][:gene_end-3] + "111" + \
						label_dict["stop_codon"][gene_end:]
				#print(label_dict["stop_codon"][gene_start-1:gene_end+1])
			else:
				pass 
				#print("Error: No STOP codon at the end of the gene, not labeling;  Investigate!!!\n", row) #TODO:  add file name to error message

			continue
			# assign the upstream and downstream regions and intragenic regions;
			# for now set upstream to 10 bases and downstream to 100 bases
			# Test the distance between prev_gene_stop and gene_start, to make sure that our
			# upstream and downstream regions do not overlap each other, or go into another gene's body;
			# in those cases, maybe make width of upstream and downstream regions smaller
			# will find prokaryotic polysistronic genes, that we will have to deal with
			upstream_length = 10
			downstream_length = 100
			if ((gene_start - prev_gene_stop) > 0): 
				# do not annotate intergenic, upstream, downstream for overlapping genes, 
				# but still keep their respective start & stop
				# TODO: will need to think further on what to do with overlapping/largeley-overlapping genes
		    
				# assign the intragenic label to 1 for everything between two genes
				label_dict["intragenic"] = label_dict["intragenic"][:prev_gene_stop] + \
						"1" * (gene_start-prev_gene_stop) + \
						label_dict["intragenic"][gene_start:]
		    
				if ((gene_start - prev_gene_stop) >= (upstream_length + downstream_length)):
					# set the upstream region for this gene
					label_dict["upstream"] = label_dict["upstream"][:gene_start - upstream_length] + \
							"1" * upstream_length + \
							label_dict["upstream"][gene_start:]
                			# set the downstream region for the previous gene
                			if (prev_gene_stop != 0): # don't label imaginary downstream gene
                    				label_dict["downstream"] = label_dict["downstream"][:prev_gene_stop] + \
                                        			"1" * downstream_length + \
                                        			label_dict["downstream"][prev_gene_stop + downstream_length:]
        
        		# update some variables
        		prev_gene_stop = gene_end
        		prev_coding_start_end = (0, 0) # used to set the intron regions and sanity checks
		
		continue
		if (row[2] in ("CDS","exon")):
        
        		if (not consider_coding_regions):
            			#print("skipping because not processing the gene: ", row)
            			continue
    
			coding_start = int(row[3])-1 # subtract 1 because annotation uses 1-based coordinates
			coding_end = int(row[4])
			coding_length = coding_end - coding_start
	
			if (coding_end > chrom_length):
				continue

			# sometimes there is a CDS and an identical exon with the same coordinates
			# sometimes there is just a CDS, sometimes there is just an exon
			# for now, skip an exon or CDS if it is identical to the previous one
			if (coding_start == prev_coding_start_end[0] and coding_end == prev_coding_start_end[1]):
	   			 #print("processing same exon/CDS as previous ", row)
			    continue
        
      			# now check for overlapping coding regions within a gene;
   			# this should not happen much in prokaryotes;  when we process eukaryotes with
        		# alternative splicing, we may need to consider this, and process things a transcript at a time
        		# allowing overlaps;  for now, lets not consider this
        		# sometimes, it can be exon1, exon2, CDS1, CDS2, which causes this problem - can get rid of this
        		# error message with more checks in that case
        		if (prev_coding_start_end[1] > coding_start):
            			#print("Diagnostic message: coding region overlapping with previous one, will ignore; investigate: ", row) #TODO: add file name
            			continue  
        
        		# assign the exon label to 1's for the stretch of this exon/coding region
        		label_dict["exon"] = label_dict["exon"][:coding_start] + "1" * coding_length \
                        	    + label_dict["exon"][coding_start+coding_length:]
            
        		# assign the intron label to 1 for everything between two exon/coding regions
        		# this will not assign the intron region after the last exon;  probably that makes sense 
        		if (prev_coding_start_end[1] != 0): # this is the first exon of the gene
            			label_dict["intron"] = label_dict["intron"][:prev_coding_start_end[1]] + \
                        		"1" * (coding_start-prev_coding_start_end[1]) + \
                                        label_dict["intron"][coding_start:]
        
        		# label all the amino acids
        		# first verify that the coding region is %3
        		# if it is not, we have some exons that start out of frame, so we will need extra calculations
        		# to figure out the "frame" for each exon
        		if (coding_length%3 != 0):
            			print("Coding region is not a multiple of 3, will not assign amino acids for now!!!: ", row, "\n")
        		else:
            			# for the exon/coding region, assign a 
            			for i in range(coding_start, coding_end, 3):
                			codon = chrom_letters[current_chrom][i:i+3]
                    			label_dict["frame1"] = label_dict["frame1"][:i] + "100" + label_dict["frame1"][i+3:]
                			if (codon in codon_to_aa_dict): 
                    				# if the codon does not map to a coding AA, such as a stop codon, or anything with N's
                    				# this will simply not assign a label
                    				aa = codon_to_aa_dict[codon]
                    				label_dict[aa] = label_dict[aa][:i] + "1" * 3 + label_dict[aa][i+3:]
                    				#print(codon_to_aa_dict[codon])
           			#print(label_dict["M"][coding_start-3:coding_end+3])
                        
        		# update some variables
        		prev_coding_start_end = (coding_start, coding_end)


	"""
	# here, output to file all the labels for this region
        # TODO:  output just the sequence, or perhaps do the three-mer or hepta-mer conversion here???
	# TODO: caveat - will not print out the "last" region with this loop; repeat under the loop?
	file_write.write(">" + current_chrom + "\n")
	file_write.write(chrom_letters + "\n")
	# here is is important that we keep the order of printing out the label_dict consistent
	# this will be the order in which we will interpret the labels in the verification pipeline
	file_write.write(label_dict["gene"] + "\n")
	file_write.write(label_dict["start_codon"] + "\n")
	file_write.write(label_dict["stop_codon"] + "\n")
	file_write.write(label_dict["upstream"] + "\n")
	file_write.write(label_dict["downstream"] + "\n")
	file_write.write(label_dict["intragenic"] + "\n")
	file_write.write(label_dict["exon"] + "\n")
	file_write.write(label_dict["intron"] + "\n")
	for aa in aa_list: # in sorted alphabetical order
		file_write.write(label_dict[aa] + "\n")
	"""
	#label_dict["intragenic"] = label_dict["intragenic"][:prev_gene_stop] + \
	#		"1" * (len(chrom_letters[current_chrom])-prev_gene_stop)
	assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))
	region_dict[current_chrom] = label_dict

	print "Finish Annotating Region: ", (current_chrom, len(label_dict["gene"]), len(chrom_letters[current_chrom]), len(region_dict[current_chrom]["sequence"]))
	return region_dict

				

for i in range(begind, min(endind, len(gff_file_list))):
	file_print_count = 0

	gff_file_name = gff_dir_name + '/' + gff_file_list[i]
	gff_file_read = gzip.open(gff_file_name, 'rb')
	#with gzip.open(gff_file_name, 'rb') as zip_file:
	#	gff_file_content = zip_file.read()
	#	#print file_content
	fna_file_name = fna_dir_name + '/' + gff_file_list[i][0 : len(gff_file_list[i]) - 6] + 'fna.gz'
	with gzip.open(fna_file_name, 'rb') as zip_file:
		fna_file_content = zip_file.read()
		#print file_content

	#gff_lines = gff_file_content.split('\n')
	fna_lines = fna_file_content.split('\n')
	#if debug == 99:
	print 'Reading GFF File: ', gff_file_name
	print 'Reading FNA File: ', fna_file_name
	#print 'Length of File: ', (len(gff_file_content), len(fna_file_content))
	#print 'Total GFF Lines: ', len(gff_lines)
	print 'Total FNA Lines: ', len(fna_lines)

	chromosome_toseq = {}
	#annotation_toseq = {}
	gene_tags = {}
	sequence_string = ""
	sequence_name = "" 

	fna_lines.append(">>>")
	for l in range(0, len(fna_lines)):
		line = fna_lines[l]
		if len(line) == 0:
			continue

		if line[0] == ">":
			if len(sequence_name) > 0 and len(sequence_string) > 0:
				chromosome_toseq[sequence_name] = re.sub(r'([^NACGT])', "N", sequence_string.upper())
				#annotation_string = ""
				#for k in range(0, len(sequence_string)):
				#	annotation_string += "0"
				#annotation_toseq[sequence_name] = annotation_string
				gene_tags[sequence_name] = []
			#print "Sequence Name: ", line, line.split(" ")[0]
			sequence_name = line.split(" ")[0][1:]
			sequence_string = ""
		else:
			sequence_string += line.strip()

	#assert(len(chromosome_toseq) == len(annotation_toseq))
	#send-> chrome_toseq, gff_file_read
	#try:
	region_dict = annotate_sequence_following_gff(chromosome_toseq, gff_file_read)
	#except:
	#	print "Something Wrong Happened"
        #IGNORING STRAND +/-
	if debug == 99:
        	print "Writing File: ", (i, gff_file_list[i])

	for key in chromosome_toseq.keys():
		if region_dict.has_key(key) == False:
			continue 

		if file_print_count % batch_size == 0 or file_write == None:
			if not file_write == None:
				file_write.close()
			file_write = gzip.open(output_dir_name + "/batch_" + str(i) + "_" + str(file_print_count) + ".gz", "w")
		
                reference_name = key
		label_dict = region_dict[reference_name]
		sequence_string = label_dict["sequence"]
		print "Label Gene: ", (reference_name, len(label_dict["gene"]), len(region_dict[key]["sequence"]), len(sequence_string))
		if not (len(label_dict["gene"]) == len(sequence_string)):
			assert(False)

		file_write.write('>' + reference_name + "," + gff_file_list[i][0 : len(gff_file_list[i]) - 7] + '\n')
                file_write.write(sequence_string + '\n')

		assert(len(label_dict["gene"]) == len(sequence_string))
		file_write.write(label_dict["gene"] + "\n")
		assert(len(label_dict["start_codon"]) == len(sequence_string))
		file_write.write(label_dict["start_codon"] + "\n")
		assert(len(label_dict["stop_codon"]) == len(sequence_string))
		file_write.write(label_dict["stop_codon"] + "\n")
		"""
		assert(len(label_dict["upstream"]) == len(sequence_string))
		file_write.write(label_dict["upstream"] + "\n")
		assert(len(label_dict["downstream"]) == len(sequence_string))
		file_write.write(label_dict["downstream"] + "\n")
		assert(len(label_dict["intragenic"]) == len(sequence_string))
		file_write.write(label_dict["intragenic"] + "\n")
		assert(len(label_dict["exon"]) == len(sequence_string))
		file_write.write(label_dict["exon"] + "\n")
		assert(len(label_dict["intron"]) == len(sequence_string))
		file_write.write(label_dict["intron"] + "\n")
		assert(len(label_dict["frame1"]) == len(sequence_string))
		file_write.write(label_dict["frame1"] + "\n")
		for aa in aa_list: # in sorted alphabetical order
			assert(len(label_dict[aa]) == len(sequence_string))
			file_write.write(label_dict[aa] + "\n")
		"""
		"""
		for label in label_dict.keys():
			annotation_string = label_dict[label]
			assert(len(annotation_string) == len(sequence_string))
			file_write.write(annotation_string + '\n')
		"""
		file_print_count += 1
	
	#if flag == True:
	#break

if not file_write == None:
	file_write.close()

print 'Total Gene found: ', count
print "Wrong index found: ", wrong


