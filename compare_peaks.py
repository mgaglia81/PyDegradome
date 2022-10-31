## this file can be used to compare two output files from PyDegradome
## It assumes that the name of the input files is TestvControl_parameters.txt
## as standard out it outputs a table of shared peaks and as standard error it outputs numbers of total peaks detected, shared peaks and peaks unique to one sample

from sys import *

if len(argv) < 7:
    print("Need arguments:\n1 = suffix test 1\n2 = suffix ctl 1\n3 = suffix test 2\n4 = suffix ctl 2\n5 = other file suffix\n6 = distance\nSyntax: python compare_peaks.py <test suffix 1> <control suffix 1> <test suffix 2> <control suffix 2> <other file suffix> <maximum distance> > <output file name>")
    exit()

## set up a class to store the data
class peak:
    def __init__(self, c):
        self.c = c
        self.peaks = {}

## set up dictionary for that class
pd = {}

## set up list of files
f_l = []

suff = [argv[1],argv[3]]
suff_comp = [argv[2],argv[4]]

for x in [0,1]:
    fn = suff[x]+"v"+suff_comp[x]+"_"+argv[5]+".txt"
    f_l.append(fn)

## set up max distance between the same peaks
diff = int(argv[6])

## go through files and record peaks

chrm_list = [] # record all chrm identifier to later use in output
ex_l = ["[","#"] # lines to exclude start with these
tot = [0,0] # number of peaks in each file

for x in [0,1]:
    stderr.write(f_l[x]+"\n")
    f = open(f_l[x], "r")
    for l in f:
        if l[0] not in ex_l:
            tot[x] += 1
            l = l.strip()
            sl = l.split("\t")
            c = sl[0]
            if c not in chrm_list:
                chrm_list.append(c)
            st = sl[1]
            c_id = c+"*"+st
            if c_id not in pd:
                pd[c_id] = peak(c_id)
            p = int(sl[4]) ## position with max count within peak = putative cut site
            h = int(sl[7]) ## max count within peak - need to make sure reading in right column
            if x == 0:  ## when going through first file record all peaks
                pd[c_id].peaks[p]={}
                pd[c_id].peaks[p][suff[x]] =[p,h]
            elif x == 1: ## when going through second file have to do more
                p_l = []
                for y in range(-diff, +diff+1): ## add potential distance to generate list of peak positions that would be considered a match
                    p_l.append(p+y)
                flag = 0
                for p_add in p_l:
                    if p_add in pd[c_id].peaks: ## check whether any of these potential match positions are already in the dictionary (i.e. were in file 1)
                        pd[c_id].peaks[p_add][suff[x]]= [p,h]
                        flag = 1
                        break
                if flag == 0: ## if no match add a new peak
                    pd[c_id].peaks[p]={}
                    pd[c_id].peaks[p][suff[x]]= [p,h]

## output

for x in [0,1]:
    mess1 = "sample"+suff[x]+": "+str(tot[x])+" peaks\n"
    stderr.write(mess1)

tot_c = 0
both = 0
firstonly = 0
secondonly = 0

chrm_list.sort()

st_list = ["+","-"]

print("#chr\tstrand\tmax_peak_"+suff[0]+"\tmax_peak_"+suff[1]+"\tcount_at_peak_"+suff[0]+"\tcount_at_peak_"+suff[1])

for c in chrm_list:
    for st in st_list:
        c_id = c+"*"+st
        if c_id in pd:
            order = list(pd[c_id].peaks.keys())
            order.sort()
            for p in order:
                tot_c += 1
                if len(pd[c_id].peaks[p]) ==2:
                    both += 1
                    new_line = [c,st,str(pd[c_id].peaks[p][suff[0]][0]),str(pd[c_id].peaks[p][suff[1]][0]),str(pd[c_id].peaks[p][suff[0]][1]),str(pd[c_id].peaks[p][suff[1]][1])]
                    print("\t".join(new_line))
                else:
                    if suff[0] in pd[c_id].peaks[p]:
                        firstonly +=1
                    elif suff[1] in pd[c_id].peaks[p]:
                        secondonly +=1

stderr.write("tot diff peaks: "+str(tot_c)+"\n")
stderr.write("peaks in both: "+str(both)+"\n")
stderr.write("peaks in "+suff[0]+" only: "+str(firstonly)+"\n")
stderr.write("peaks in "+suff[1]+" only: "+str(secondonly)+"\n")
