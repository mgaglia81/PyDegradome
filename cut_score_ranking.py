from sys import *
from Bio import SeqIO
import matplotlib.pyplot as plt
import numpy as np
import re

if len(argv) < 7 :
    print "Need arguments:\n1 = peak file\n2 = gtf file\n3 = motif pwm\n4 = genome fasta\n5 = ouput pdf name\n6 = position of nt 0 relative to total tested sequence length\nSyntax: cut_score_ranking.py <peak file> <gtf file> <pwm file> <genome fasta> <output pdf name> <pos 0>"
    exit()

## define motif loading and scoring functions

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
    for x in range(0, len(R)-len_m):
        s = R[x:x+len_m]
        p = 0 ## use p = 1 if probabilities
        if "N" in s:
            continue
        for y in range(0,len(s)):
            nuc = s[y]
            if nuc not in pwm:
                nuc = conv_d[nuc]
            val = pwm[nuc][y]
            p += val ## use p = p*val if probabilities
        score_l.append(p) ## use math.log(p,2 or 10) if using probabilities
    return score_l

## load peaks into a dictionary

pd = {}

pf = open(argv[1],"r")

tot_p = 0

class peak:
    def __init__(self, p, l):
        self.p = p
        self.l = l
        self.t_id = []
        self.ex = {}

for l in pf:
    if l[0] != "#" and l[0] !="[":
        tot_p += 1
        l = l.strip()
        sl = l.split("\t")
        pl = sl[0:6]
        gr = sl[0]+"*"+sl[1]
        if gr not in pd:
            pd[gr] = {}
        p_id = int(sl[2]) ## appropriate if using only cut sites that have exact same position in both datasets
        pd[gr][p_id] = peak(p_id,pl)

## read in gtf and identify:transcript it comes from, all exons belonging to transcript

gtf = open(argv[2],"r")

class transcript:
    def __init__(self, t_id, strand):
        self.t_id = t_id
        self.exlist = []
        self.strand = strand
        self.seq = ""

td = {}

for l in gtf:
    l = l.strip()
    sl = l.split("\t")
    gr = sl[0]+"*"+sl[6]
    strand = sl[6]
    if gr in pd:
        if gr not in td:
            td[gr] = {}
        if sl[2] == "exon":
            t_id_s = re.search("; transcript_id \"(.*?)\";", sl[8])
            t_id = t_id_s.group(1)
            if t_id not in td[gr]:
                td[gr][t_id] = transcript(t_id, strand)
            td[gr][t_id].exlist.append([int(sl[3]), int(sl[4])])
            for p in pd[gr]:
                if int(sl[3])<=p <=int(sl[4]):
                    pd[gr][p].ex[t_id] =[int(sl[3]), int(sl[4])]
                    pd[gr][p].t_id.append(t_id)

## load motif

m = open(argv[3],"r")

pwm = load_motif(m)
stderr.write("motif loaded\n")
len_seq = len(pwm["A"])
stderr.write("length of sequence: "+str(len_seq)+"\n")
align_num = int(argv[6])

## read in genome sequence

g_seq_d = {}

fs = open(argv[4], "rU")

for record in SeqIO.parse(fs,"fasta"):
    a = record.id
    b = record.seq
    g_seq_d[a] = b

## recreate transcripts based on annotation, find position of p to compare, score sequences

tot = 0
yes = 0
attopten = 0
attoptwenty = 0
rel_match_l = []

for gr in pd:
    temp2 = gr.split("*")
    chrm = temp2[0]
    strand = temp2[1]
    for p in pd[gr]:
        tot += 1
        t_l = pd[gr][p].t_id
        rel_match = 2.0
        close = "NA"
        for z in range(0, len(t_l)):
            t_id = t_l[z]
            ex_p = pd[gr][p].ex[t_id]
            exlist = td[gr][t_id].exlist
            seqRNA = ""
            pos = 0
            for ex in exlist:
                beg = ex[0]-1
                end = ex[1]
                seq = g_seq_d[chrm][beg:end]
                if strand == "+":
                    seqRNA += seq
                    if ex == ex_p:
                        temp3 = p - beg - 1
                    elif beg < ex_p[0]:
                        temp3 = end - beg
                    else:
                        temp3 = 0
                elif strand == "-":
                    seq = seq.reverse_complement()
                    seqRNA = seq+seqRNA
                    if ex == ex_p:
                        temp3 = end - p
                    elif beg > ex_p[0]:
                        temp3 = end - beg
                    else:
                        temp3 = 0
                pos += temp3
            score = m_seq_score(pwm, seqRNA, len_seq)
            score_ord = [x for x in score]
            score_ord.sort(reverse=True)
            ord_pos_ind = score_ord.index(score[pos-align_num]) + 1
            rel_pos = float(ord_pos_ind)/float(len(score))
            if rel_pos < rel_match:
                rel_match = rel_pos
            top_score_twenty = []
            top_score_ten = []
            for y in range(0,21):
                val = score_ord[y]
                ind = score.index(val)
                top_score_twenty.append(ind)
                if y < 11:
                    top_score_ten.append(ind)
            pos_max = score.index(max(score))
            extra = []
            for w in range(pos-align_num-2,pos-align_num+3):
                extra.append(str(score[w]))
            extra_str = ",".join(extra)
            #new_line = [chrm, strand, str(p), str(pos), str(seqRNA[pos-10:pos+10]), str(pos_max), str(max(score)), extra_str]
            #print "\t".join(new_line) ## option to print out to standard ouput
            if pos == pos_max+align_num:
                if close != "yes":
                    close = "yes"
            elif pos-align_num in top_score_ten:
                if close != "yes":
                    close = "attopten"
            elif pos-align_num in top_score_twenty:
                if close != "yes" and close != "attopten":
                    close = "attoptwenty"
        if close == "yes":
            yes += 1
        elif close == "attopten":
            attopten += 1
        elif close == "attoptwenty":
            attoptwenty += 1
        rel_match_l.append(rel_match)

plt.figure(1)
plt.hist(rel_match_l,cumulative = True, histtype = "step", normed = True, color = "b", lw = 3)
plt.savefig(argv[5], format="PDF")

message = "\ntot RNAs = %d\nmax score at cut site = %d\ncut site score within top 10 scores = %d\ncut site score within top 20 scores = %d\n" % (tot, yes, attopten, attoptwenty)
stderr.write(message)
