# PyDegradome, a program to find cuts site in degradome-seq/PARE experiments
# with control and test datasets
#
# Authors   : Marta Gaglia (UW-Madison), Chris H. Rycroft (UW-Madison/LBL)
# Email     : marta.gaglia@wisc.edu, rycroft@wisc.edu

from sys import *
from math import *
from bisect import bisect
from scipy import interpolate
import numpy as np
from thr_class import *
import re

# Default value for the total number of sites that could possibly be cut, based
# on 1% of the human genome
tot_sites = 32348300.

# Error message to give if command-line arguments can't be correctly parsed
def error_message():
    stderr.write("PyDegradome: Unrecognized command-line options; type \"python PyDegradome.py -h\"\nfor more information.\n")
    exit()

# Help message
def help_message():
    print("Syntax:\n" \
          "python PyDegradome.py [options] <gtf file> <control .sam file> <test .sam file>\n" \
          "       <confidence level> <window size> <multiplicative factor> <output file name>\n\n" \
          "Options:\n" \
          "  -h, --help        Show the syntax of the program and exit\n" \
          "  -t <total sites>  Override the default value of total possible read sites")
    exit()

# Check if help requested
if len(argv) > 1 and (argv[1] == "-h" or argv[1] == "--help"):
    help_message()

# Check for enough command-line arguments
if len(argv) < 8:
    error_message()

# Parse any options
while len(argv) > 8:
    cmd=argv.pop(1)
    if cmd == "-h" or cmd == "--help":
        help_message()
    elif cmd == "-t":
        if argv == 8:
            error_message()
        tot_sites=float(argv.pop(1))
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

        # key = gene indentifier, value = list of hit windows [start,stop];
        # position in list: 0 = + strand, 1 = - strand
        self.hit = [{},{}]

        # key = gene indentifier, value = list of hit windows [start,stop];
        # position in list: 0 = + strand, 1 = - strand
        self.hitrefined = [{},{}]
        self.order = [[],[]]

    # Collapses alternative transcripts into a big one
    def add_genes(self,strand, start, stop):
        if strand == "+":
            ind = 0
        elif strand == "-":
            ind = 1
        if len(self.start[ind]) > 0:
            if start < self.start[ind][-1]:
                print(start, self.start[ind][-1])
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
    def add_hit(self,s,ex_id,p1,p2):
        if ex_id not in self.hit[s]:
            self.hit[s][ex_id]=[]
        self.hit[s][ex_id].append([p1,p2])

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
# First generate a list where pos = read count, value = freq of that read count
# 0 = control, 1 = test
chd = [{},{}]

# Scanning window size to integrate over
w = int(argv[5])

# Record the number of windows with specific read count as a dictionary with
# key = read count, value = # of instances (one for test, and one for control
# samples)
for chrm in gd:
    for x in [0,1]:
        for ex in gd[chrm].gene_l[x]:
            for z in [0,1]:
                temp_hl = gd[chrm].gene_l[x][ex][z]
                for nt in range(0, len(temp_hl)-w+1):
                    w_hl = sum(temp_hl[nt:nt+w])
                    ch = l_scale(w_hl)
                    if ch != 0:
                        if ch not in chd[z]:
                            chd[z][ch] = 0
                        chd[z][ch] += 1

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
    fn = argv[7] + "_" + fn_l[x] + "_heightfreq.txt"
    fh = open(fn, "w")
    for y in range(len(chd_l[x])):
        nl = [str(y),str(chd_l[x][y])]
        nlj = "\t".join(nl)
        fh.write(nlj+"\n")
    fh.close()

# Generate the threshold table using the lists above with the given input
# confidence level (stored in argv[4])
th=thresh_tab()
th.fit_prior(chd_l,tot_sites)
th.make_table(float(argv[4]))
th.save_table(argv[7]+"_test.tab")

# Go through counts in genes, identify count in controls, test whether count in
# test exceed the threshold by looking up the threshold in the generated table,
# and record windows that passed
count = 0
co = float(argv[6])

stderr.write("finding peaks:\n")

for chrm in gd:
    for x in [0,1]:
        for ex in gd[chrm].gene_l[x]:
            count += 1
            if count % 100000 == 0:
                stderr.write("processed "+str(count)+" genes\n")
            flag = 0
            for ds in [0,1]:
                if gd[chrm].genecov[x][ex][ds] > 10:
                    flag += 1
            if flag == 2:
                ge_cov_ratio = float(gd[chrm].genecov[x][ex][1])/float(gd[chrm].genecov[x][ex][0])
                factor = co*ge_cov_ratio

                # If the factor is larger 100, it will be artificially set to
                # 100 and a warning message will be printed
                if factor > 100:
                    err_mess = "Ratio of " + str(factor) + " found, chrm "+ chrm+", gene "+ex
                    stderr.write(err_mess+"\n")
                    factor = 100
                temp_t = gd[chrm].gene_l[x][ex][1]
                temp_c = gd[chrm].gene_l[x][ex][0]
                for nt in range(0, len(temp_c) - w + 1):
                    w_t = sum(temp_t[nt:nt+w])
                    w_c = sum(temp_c[nt:nt+w])
                    if w_t != 0 or w_c != 0:
                        stat = th.thr(w_c,factor)
                        if w_t >= stat and max(temp_t[nt:nt+w]) > 1:
                            gd[chrm].add_hit(x,ex,nt,nt+w)

# Trim windows so that 0 or 1 count positions are eliminated. For example,
# [0,1,2,4,5] becomes [2,4,5].
for c in gd:
    for x in [0,1]:
        for ex in gd[c].hit[x]:
            temp = gd[c].hit[x][ex]
            cov = gd[c].gene_l[x][ex][1]
            if len(temp) != 0:
                for i in range(0,len(temp)):
                    a = cov[temp[i][0]:temp[i][1]]
                    c1 = 0
                    c2 = 0
                    for j in range(0,len(a)):
                        if a[j] < 2:
                            c1 += 1
                        else:
                            break
                    for j in range(len(a) - 1,-1,-1):
                        if a[j] < 2:
                            c2 += 1
                        else:
                            break
                    gd[c].hit[x][ex][i]=[temp[i][0]+c1,temp[i][1]-c2]

# Merge windows into bigger peaks
for c in gd:
    for x in [0,1]:
        for ex in gd[c].hit[x]:
            temp_hit = gd[c].hit[x][ex]
            if len(temp_hit) != 0:
                prev_h = temp_hit[0][0]
                temp_start = [prev_h]
                temp_end = []
                for i in range(1,len(temp_hit)):
                    next_h = temp_hit[i][0]
                    if next_h > prev_h+10:
                        temp_start.append(next_h)
                        temp_end.append(temp_hit[i-1][1])
                    prev_h = next_h
                temp_end.append(temp_hit[-1][1])
                for j in range(0,len(temp_start)):
                    if ex not in gd[c].hitrefined[x]:
                        gd[c].hitrefined[x][ex] = []
                    gd[c].hitrefined[x][ex].append([temp_start[j],temp_end[j]])

# Construct main output file
#
# Format is 3 lines per peak:
#
# The first line includes chr_gene-start_gene-end, strand, position of peak start, position of peak
# end, the position with the highest read count, the ratio between all reads in
# the genes in test vs. control, the factor used for thresholding (user-defined
# ratio * ratio of gene reads), the read count at the position listed as
# "position with the highest read count" in test and control sample, the total
# read count in the peak in the test and control samples.
#
# The second and third lines print out the list of counts for positions
# starting 10 nt 5' of peak start and ending 10 nt 3' of peak end (fewer if the
# peak is less than 10 nt away from the end of an annotated gene) for each of
# control and test dataset.

# Open file and print header line
outf = open(argv[7],"w")
outf.write("#chr_gene-start_gene-end\tstrand\tpeak_start\tpeak_stop\tmax_peak_test\tratio\tfactor\tmax_count_test\tmax_count_ctl\ttot_count_test\ttot_count_ctl\n")

chrm_list = list(gd.keys())
chrm_list.sort()

strand_l = ["+","-"]

for c in chrm_list:
    for x in [0,1]:
        order = list(gd[c].hitrefined[x].keys())
        order.sort()
        for ex in order:
            temp_fin = gd[c].hitrefined[x][ex]
            cov_t = gd[c].gene_l[x][ex][1]
            cov_c = gd[c].gene_l[x][ex][0]
            info = ex.split("-")
            ex_start = int(info[0])
            ex_end = int(info[1])
            strand = strand_l[x]
            if len(temp_fin) != 0:
                for h in temp_fin:
                    test = cov_t[h[0]:h[1]]
                    control = cov_c[h[0]:h[1]]
                    maxtest = max(test)
                    maxcontrol = max(control)
                    peak_rel = np.where(test == maxtest)[0][0]
                    start = ex_start + h[0]
                    end = ex_start + h[1]
                    peak = start + peak_rel
                    tot_t = sum(test)
                    tot_c = sum(control)
                    ratio = float(gd[c].genecov[x][ex][1])/float(gd[c].genecov[x][ex][0])
                    factor = ratio*co
                    chrm_id = c + "_" + str(ex_start) + "_" + str(ex_end)
                    new_line = [chrm_id, strand, str(start), str(end), str(peak), str(ratio), str(factor), str(maxtest), str(maxcontrol),str(tot_t),str(tot_c)]
                    nlj2 = "\t".join(new_line)
                    outf.write(nlj2+"\n")
                    if h[0]-10 > 0 :
                        outf.write(str(cov_c[h[0]-10:h[1]+10])+"\n")
                        outf.write(str(cov_t[h[0]-10:h[1]+10])+"\n")
                    else:
                        outf.write(str(cov_c[0:h[1]+10])+"\n")
                        outf.write(str(cov_t[0:h[1]+10])+"\n")

outf.close()
