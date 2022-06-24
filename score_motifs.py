### this program can be used to score motifs in a fasta file against a position weight matrix (PWM); takes one or more fasta files and assumes length of sequences in fasta match the size of the PWM

from sys import *
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

if len(argv) < 4:
    print "Need arguments:\n1 = pwm\n2 = title of pdf output\n3 and following = fasta to score\nSyntax: python score_motifs.py <pwm> <title of pdf output file> <fasta to score 1> <fasta to score 2> ..."
    exit()

## function to load in fasta file - used the identifier after ">" as id and the string on the following line as the sequence
def load_fasta(f):
    d = {}
    for l in f:
        l = l.strip()
        if l[0] == ">":
            t_id = l[1:]
        else:
            d[t_id]= l
    return d

## function to load in motif - motif needs to be in following format: rows = nucleotide, column = position and be in log likelihood - input is a text file, first row can be headers if "#" in front
def load_motif(f):
    lm = {}
    for l in f:
        if l[0] != "#":
            l = l.strip()
            sl = l.split("\t")
            float_sl = [float(x) for x in sl[1:]]
            lm[sl[0]]=float_sl
    return lm

## function to score motif - takes in dictionary where values are sequences, and motif loaded with load_motif() function - returns a list of the scores for each sequence
def m_seq_score(pwm, d):
    score_l = []
    for k in d:
        s = d[k]
        p = 0 # would do p = 1 if using likelihood instead
        if "N" in s:
            continue
        for y in range(0,len(s)):
            nuc = s[y]
            val = pwm[nuc][y]
            p = p+val ## would do p = p*val if using likelihood instead
        score_l.append(p) ## would use math.log(p,10) if using likelihood instead
    return score_l

## load motif file as first argument
m = open(argv[1],"r")

pwm = load_motif(m)
print "pwm loaded"

## prepare settings for plotting a histogram - manually enter values to use for colors and legends
plt.figure(1)
colors = ["DarkOrchid","g","y","b"]
labels = ["sample1","sample2", "sample3","sample4"]

all_scores =[]

print "\nanalyzing files:"

for x in range(3,len(argv)): ## takes multiple files if needed
    print argv[x] ## prints out name of each file as confirmation
    fs = open(argv[x],"r")
    ds = load_fasta(fs)
    lscore = m_seq_score(pwm,ds)
    all_scores.append(lscore)
    #plt.axis([-15,15,0,0.25])
    plt.hist(lscore, normed = True, color = colors[x-3], label = labels[x-3], alpha = 0.5)

plt.legend(loc=0)
plt.savefig(argv[2], format="PDF")

print "\nKS test p-values"

if (len(argv) - 3) > 1:
    for x in range(len(all_scores)-1):
        for y in range(x+1,len(all_scores)):
            print argv[x+3],"vs.",argv[y+3]
            a = all_scores[x]
            b = all_scores[y]
            KS_stats = stats.ks_2samp(a,b)
            pval = KS_stats[1]
            print pval
