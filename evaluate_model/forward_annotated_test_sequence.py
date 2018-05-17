from os import walk
import gzip
import csv
import re

#comm = MPI.COMM_WORLD
#name=MPI.Get_processor_name()
#print("hello world")
#print(("name :",name,"my rank is: ",comm.rank,"out of total: ",comm.size))

genmark_organisms = ['Bacillus subtilis', 'Escherichia coli', 'Haemophilus influenzae', \
            'Helicobacter pylori', 'Mycoplasma genitalium', 'Mycoplasma pneumoniae', \
            'Methanobacterium thermoauthotrophicum', 'Synechocystis', \
            'Archeoglobus fulgidus', 'Marinobacterium jannaschii']

file_read = open("../dataset/assembly_summary.txt", "r")
complete = {}
unique_species = {}
for line in file_read:
    columns = line.split("\t")
    if len(columns) < 21:
        print(line)
        continue

    #species_prefix = columns[7].split(" ")[0]
    genmark_flag = False
    for k in range(0, len(genmark_organisms)):
        if columns[7].startswith(genmark_organisms[k]):
            genmark_flag = True
            species_prefix = genmark_organisms[k]
            break
    if genmark_flag == False:
        continue

    if species_prefix in unique_species:
        #print(columns)
        continue
    else:
        if columns[11] == "Complete Genome":
            #print columns
            print(columns[0], columns[7], columns[11])
            complete[columns[0]] = True
            unique_species[species_prefix] = True
print("Length of Complete Genome = ", len(complete))


fna_dir_name = "../dataset/GbBac_FNA"#"debug_files/"#"dataset/GbBac_FNA"#'test_fna'
gff_dir_name = "../dataset/GbBac_GFF"#"debug_files/"#"dataset/GbBac_GFF"#'test_gff'
output_dir_name = "./"#"./batch_files"#"./test_files"

gff_file_list = []
for (dirpath, dirnames, filenames) in walk(gff_dir_name):
    for i in range(0, len(filenames)):
        if ("GCA_" + filenames[i].split("_")[1]) in complete and "gff" in filenames[i]:
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
            if (current_chrom != ""): 
                # ****we are done with the whole region - output everything
                # here, output to file all the labels for this region
                # TODO:  output just the sequence, or perhaps do the three-mer or hepta-mer conversion here???
                # TODO: caveat - will not print out the "last" region with this loop; repeat under the loop?

                assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))
                region_dict[current_chrom] = label_dict
                print("Finish Annotating Region: ", (current_chrom, len(label_dict["gene"]), \
                        len(chrom_letters[current_chrom]), len(region_dict[current_chrom]["sequence"])))
        
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
            
            assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))

        if (row[2] == "gene"):
        
            # Important: for now, only consider protein-coding genes 
            # - tRNA genes etc do not start with ATG, or end with stop
            # TODO:  do sanity check here - should be keeping the majority of genes
            if (("gene_biotype=protein_coding" not in row[8]) or (row[6] == "-")):
                #print("skipping: ", row)
                continue
        
            #print(row)
            gene_start = int(row[3])-1 # subtract 1 because annotation uses 1-based coordinates
            gene_end = int(row[4])
            gene_length = gene_end - gene_start

            #if (gene_start > 30000):#Condition for debugging
            #    break

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

            # check that the stop codon is the end of the gene and set the label
            if (chrom_letters[current_chrom][gene_end-3:gene_end] in ("TAG", "TAA", "TGA")) :
                label_dict["stop_codon"] = label_dict["stop_codon"][:gene_end-3] + "111" + \
                                            label_dict["stop_codon"][gene_end:]
                #print(label_dict["stop_codon"][gene_start-1:gene_end+1])

    
    assert(len(label_dict["gene"]) == len(chrom_letters[current_chrom]))
    region_dict[current_chrom] = label_dict

    print("Finish Annotating Region: ", (current_chrom, len(label_dict["gene"]), len(chrom_letters[current_chrom]), len(region_dict[current_chrom]["sequence"])))
    return region_dict

                

for i in range(begind, min(endind, len(gff_file_list))):
    file_print_count = 0

    gff_file_name = gff_dir_name + '/' + gff_file_list[i]
    gff_file_read = gzip.open(gff_file_name, 'rt')
    #with gzip.open(gff_file_name, 'rb') as zip_file:
    #    gff_file_content = zip_file.read()
    #    #print file_content
    fna_file_name = fna_dir_name + '/' + gff_file_list[i][0 : len(gff_file_list[i]) - 6] + 'fna.gz'
    with gzip.open(fna_file_name, 'rt') as zip_file:
        fna_file_content = zip_file.read()
        #print file_content

    #gff_lines = gff_file_content.split('\n')
    fna_lines = fna_file_content.split('\n')
    #if debug == 99:
    print('Reading GFF File: ', gff_file_name)
    print('Reading FNA File: ', fna_file_name)
    #print 'Length of File: ', (len(gff_file_content), len(fna_file_content))
    #print 'Total GFF Lines: ', len(gff_lines)
    print('Total FNA Lines: ', len(fna_lines))

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
                #    annotation_string += "0"
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
    #    print "Something Wrong Happened"
        #IGNORING STRAND +/-
    if debug == 99:
        print("Writing File: ", (i, gff_file_list[i]))

    for key in chromosome_toseq.keys():
        if not key in region_dict:
            continue 

        if file_print_count % batch_size == 0 or file_write == None:
            if not file_write == None:
                file_write.close()
            file_write = gzip.open(output_dir_name + "/forward_batch_" + str(i) + "_" + str(file_print_count) + ".gz", "wt")
        
        reference_name = key
        label_dict = region_dict[reference_name]
        sequence_string = label_dict["sequence"]
        print("Label Gene: ", (reference_name, len(label_dict["gene"]), len(region_dict[key]["sequence"]), len(sequence_string)))
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
        
        file_print_count += 1
    
    #if flag == True:
    #break

if not file_write == None:
    file_write.close()

print('Total Gene found: ', count)
print("Wrong index found: ", wrong)


