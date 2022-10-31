### this program can be used to score motifs in a fasta file against a position weight matrix (PWM); takes one or more fasta files and scans through each sequence, returns the highest score for each sequence
##
from Bio import SeqIO
from sys import *
import math
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

if len(argv) < 4:
    print("Need arguments:\n1 = pwm\n2 = title of pdf output\n3 and following = fasta to score\nSyntax: python score_motifs_wholetranscr.py <pwm> <title of pdf output file> <fasta to score 1> <fasta to score 2> ...")
    exit()

def load_motif(f):
    lm = {}
    for l in f:
        if l[0] != "#":
            l = l.strip()
            sl = l.split("\t")
            float_sl = [float(x) for x in sl[1:]]
            lm[sl[0]]=float_sl
    return lm

def m_seq_score(pwm, R, len_m): ## function to score the fasta files
    score_l = []
    for x in range(0, len(R)-len_m+1):
        s = R[x:x+len_m]
        if len(R) == len_m:
            print(s)
        p = 0 ## use p =0 if probabilities
        if "N" in s:
            continue
        for y in range(0,len(s)):
            nuc = s[y]
            val = pwm[nuc][y]
            p += val ## use p = p*val if probabilities
        score_l.append(p) ## use math.log(p,2 or 10) if using probabilities
    return score_l

m = open(argv[1],"r")

pwm = load_motif(m)
print("PWM loaded")
len_seq = len(pwm["A"])
print("length of sequence:", len_seq)

plt.figure(1)
colors = ["g","DarkOrange","r","DarkOrchid"]
labels = ["sample1", "sample2","sample3"]

all_scores =[]

print("\nanalyzing files:")

for x in range(3,len(argv)):
    print(argv[x])
    fs = open(argv[x],"r")
    d = {}
    l_score = []
    flag = 0
    for record in SeqIO.parse(fs,"fasta"):
        a = record.id
        b = record.seq
        d[a] = b
    for rec in d:
        R = d[rec]
        #print len(R)
        if len(R) >= len_seq:
            score = m_seq_score(pwm,R, len_seq)
            if len(score) == 0:
                print(len(R))
                print(score)
            score_m = max(score) ## change to min(score) if want to plot the worst scores
            l_score.append(score_m) ## change to .extend(score) if want to plot all possible scores
    plt.hist(l_score , normed = True, color = colors[x-3],  label = labels[x-3], alpha = 0.5)
    all_scores.append(l_score)

plt.legend(loc=0)
plt.savefig(argv[2], format="PDF")

print("\nKS test p-values")

for x in range(len(all_scores)-1):
    for y in range(x+1,len(all_scores)):
        print(argv[x+3],"vs.", argv[y+3])
        a = all_scores[x]
        b = all_scores[y]
        KS_stats = stats.ks_2samp(a,b)
        pval = KS_stats[1]
        print(pval)
