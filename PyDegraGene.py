# PyDegraGene, a program to find genes that are cut by RNases in degradome-seq/PARE experiments
# with control and test datasets
#
# Authors   : Marta Gaglia (Univ. of Wisconsin-Madison) Chris H. Rycroft (Univ. of Wisconsin-Madison/LBL)
# Email     : marta.gaglia@wisc.edu

from sys import *
from math import *
from bisect import bisect
from scipy import interpolate
import numpy as np
from thr_class import *
import re

# Error message to give if command-line arguments can't be correctly parsed
def error_message():
    stderr.write("PyDegradome: Unrecognized command-line options; type \"python PyDegradome.py -h\"\nfor more information.\n")
    exit()

# Help message - for CHR: changed to remove window setting
def help_message():
    print "Syntax:\n" \
          "python PyDegradome.py [options] <gtf file> <control .sam file> <test .sam file>\n" \
          "       <confidence level> <multiplicative factor> <output file name>\n\n" \
          "Options:\n" \
          "  -h, --help        Show the syntax of the program and exit\n" \
          "  -t <total genes>  Override the default value of total possible genes to use"
    exit()

# Check if help requested
if len(argv) > 1 and (argv[1] == "-h" or argv[1] == "--help"):
    help_message()

# Check for enough command-line arguments
if len(argv) < 7:
    error_message()

# Parse any options
tot_genes = 0 # for CHR: to count genes from annotation file later
while len(argv) > 7:
    cmd=argv.pop(1)
    if cmd == "-h" or cmd == "--help":
        help_message()
    elif cmd == "-t":
        if argv == 7:
            error_message()
        tot_genes=float(argv.pop(1))
    else:
        error_message()

# Define class to store annotation, peaks indexed by chromosome
class chromosomes:  # Record boundaries of gene locus
    def __init__(self,ID):
        self.ID = ID
        self.start= [[],[]]         
        self.stop = [[],[]]          

        # key = gene indentifier; value = list of reads (5') that map to that gene;
        # position in list: 0 = + strand, 1 = - strand
        self.genecov = [{},{}]          

        # key = gene indentifier, value = list of lists, one for each sample,
        # each as long as the gene length;
        # position in list: 0 = + strand, 1 = - strand
        self.gene_l = [{},{}]       

        # key = gene indentifier, value = list of read #s for ctl (0) and test (1)
        # position in list: 0 = + strand, 1 = - strand
        self.hit = [{},{}]

        self.order = [[],[]]
    
    # Collapses alternative transcripts into a big one
    def add_genes(self,strand, start, stop):            
        if strand == "+":
            ind = 0
        elif strand == "-":
            ind = 1
        if len(self.start[ind]) > 0:                
            if start < self.start[ind][-1]:         
                print start, self.start[ind][-1]            
                exit()
            if start < self.stop[ind][-1]:                  
                if stop >  self.stop[ind][-1]:              
                    self.stop[ind][-1] = stop               
            else:
                self.start[ind].append(start)               
                self.stop[ind].append(stop)                 
        else:
            self.start[ind].append(start)               
            self.stop[ind].append(stop)                 
    def add_hit(self,s,ex_id, value_list):
        self.hit[s][ex_id]=value_list

# Function to find which gene the aligned read comes from
def find_gene(ex_start_list,ex_stop_list,n):            
    poss_ex = []
    a = bisect(ex_start_list,n) - 1
    if a >= 0 and n <= ex_stop_list[a]:
        poss_ex.append(ex_start_list[a])
        poss_ex.append(ex_stop_list[a])
    return poss_ex

# Adds to read count per gene and per position
def add_cov(l,d1,d2,s,p):
    k = str(l[0])+"-"+str(l[1])
    d1[k][s] += 1
    p2 = p-l[0]
    d2[k][s][p2] += 1                        

# Set up the gd dictionary based on the gtf file and record position of genes
# in each chromosome
stderr.write("processing annotation:\n")

gtf = open(argv[1], "r")
gd = {}


for l in gtf:
    l = l.strip()
    sl = l.split("\t")
    if sl[2] == "gene":            
        chrm = sl[0]
        if chrm not in gd:
            gd[chrm] = chromosomes(chrm)
        strand = sl[6]
        start = int(sl[3])
        stop = int(sl[4])
        gd[chrm].add_genes(strand, start, stop)

# Go through the read file, identify the gene the read matches to, generate a
# list for that gene if needed and count where the 5' most end of the read
# falls
stderr.write("\nprocessing control and test read files:\n")

index_list = [argv[2], argv[3]] ## sets input file list

tcr = [0,0] ## records total # of aligned reads in input
cr = [0,0] ## records total # of reads used in analysis

p=re.compile("^\d*M$")

for s in range(len(index_list)):
    filename = index_list[s]
    stderr.write(filename+"\n")
    f = open(filename, "r")
    for l in f:
        if l[0] != "@":
            l = l.strip()
            sl = l.split("\t")
            if sl[2] != "*":                        
                pos = "not set"                     
                tcr[s] += 1                         
                for fl in range(11,len(sl)):        
                    if sl[fl][0:2] == "XS":         
                        XS = sl[fl]                 

            # For simplicity, reads with insertions or deletions relative to
            # the reference are excluded
                if "I" in sl[5] or "D" in sl[5] or "S" in sl[5]:     
                    continue

                if int(sl[1]) & 16:
                
                # Removes reads that mapped to strand complementary to the
                # reference
                    if XS == "XS:A:-": 
                        x = 1
                        if (p.search(sl[5])):
                            pos = int(sl[3]) + len(sl[9]) - 1
                        else:
                            int1 = sl[5].split("M")
                            int2 = int1[1].split("N")
                            pos = int(sl[3]) + int(int1[0]) + int(int2[0]) + int(int2[1]) - 1
                else:

                # Removes reads that mapped to strand complementary to the
                # reference
                    if XS == "XS:A:+": 
                        x = 0
                        pos = int(sl[3])
                chrm = sl[2]
                if pos != "not set" and chrm in gd:     
                    poss_gene = find_gene(gd[chrm].start[x],gd[chrm].stop[x],pos)             
                    if len(poss_gene) != 0:
                        k = str(poss_gene[0])+"-"+str(poss_gene[1])                
                        if k not in gd[chrm].gene_l[x]:                                   
                            gd[chrm].genecov[x][k] = [0,0]        
                            ge_len = poss_gene[1]-poss_gene[0]+1   
                            gd[chrm].gene_l[x][k] = []                    
                            gd[chrm].gene_l[x][k].append(np.zeros(ge_len, dtype=int))      
                            gd[chrm].gene_l[x][k].append(np.zeros(ge_len, dtype=int))      
                        add_cov(poss_gene,gd[chrm].genecov[x],gd[chrm].gene_l[x],s,pos)       
                        cr[s] += 1


tcr_s = "total reads control: "+str(tcr[0])+"\ntotal reads test: "+str(tcr[1])+"\n"
cr_s = "reads used control: "+str(cr[0])+"\nreads used test: "+str(cr[1])+"\n"
stderr.write(tcr_s)
stderr.write(cr_s)
stderr.write("\nparameters of prior distributions:\n")

# Generate a table of number of positions with a certain read count to it
#
# First generate a list where key = read count, value = freq of that read count
# 0 = control, 1 = test
chd = [{},{}] 

# Record the number of windows with specific read count as a dictionary with
# key = read count, value = # of instances (one for test, and one for control
# samples)

tot_genes_temp = 0
for chrm in gd:
    for x in [0,1]:
        for ex in gd[chrm].cov[x]: 
            tot_genes_temp +=1
            for z in [0,1]:
                ch = l_scale(gd[chrm].genecov[x][ex][z])
                if ch != 0:
                    if ch not in chd[z]: 
                        chd[z][ch] = 0
                    chd[z][ch] += 1

if tot_genes == 0: #if this was not overridden in run settings
    tot_genes = tot_genes_temp

# Convert the dictionary to a list
chd_l = [[],[]]
for x in [0,1]:
    max_val = 1+max(chd[x].keys())
    chd_l[x] = [0]*max_val
    for k in chd[x]:
        chd_l[x][k]= chd[x][k]

# Print a table of the list
fn_l = ["ctl","test"]
for x in [0,1]:
    fn = argv[-1] + "_" + fn_l[x] + "_heightfreq.txt"
    fh = open(fn, "w")
    for y in range(len(chd_l[x])):
        nl = [str(y),str(chd_l[x][y])]
        nlj = "\t".join(nl)
        fh.write(nlj+"\n")
    fh.close()

# Generate the threshold table using the lists above with the given input
# confidence level (stored in argv[4])
th=thresh_tab()
th.fit_prior(chd_l,tot_genes)
th.make_table(float(argv[4]))
th.save_table(argv[-1]+"_test.tab")

# Go through counts in genes, identify count in controls, test whether count in
# test exceed the threshold by looking up the threshold in the generated table,
# and record windows that passed
count = 0
co = float(argv[5])

stderr.write("finding genes:\n")

for chrm in gd:
    for x in [0,1]:
        for ex in gd[chrm].genecov[x]:               
            count += 1
            if count % 100000 == 0:
                stderr.write("processed "+str(count)+" genes\n")
            if gd[chrm].genecov[x][ex][0] != 0 or gd[chrm].genecov[x][ex][1] != 0:
                stat = th.thr(gd[chrm].genecov[x][ex][0],factor)
                if gd[chrm].genecov[x][1] >= stat:
                    gd[chrm].add_hit(x,ex, gd[chrm].genecov[x][ex])

# Construct main output file
#
# Each line includes chr_gene-start_gene-end, strand, count for test, count for control.

# Open file and print header line
outf = open(argv[-1]+".txt","w")
outf.write("#chr_gene-start_gene-end\tstrand\tmax_count_test\tmax_count_ctl\n")

chrm_list = gd.keys()
chrm_list.sort()

strand_l = ["+","-"]

for c in chrm_list:
    for x in [0,1]:
        order = gd[c].hit[x].keys()
        order.sort()
        for ex in order:
            info = ex.split("-")
            ex_start = int(info[0])
            ex_end = int(info[1])
            strand = strand_l[x]
            chrm_id = c + "_" + str(ex_start) + "_" + str(ex_end)                          
            new_line = [chrm_id, strand, gd[c].hit[x][ex][1], gd[c].hit[x][ex][0] ]                                             
            nlj2 = "\t".join(new_line)
            outf.write(nlj2+"\n")

outf.close()
