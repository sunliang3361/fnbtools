#!/usr/bin/env python


#encoding:utf-8
"""
VarDiff.py
Copyright (c) 2017 Noble Research Institute, LLC.

Created by Liang Sun on 2016-10-13
Updated and maintained by Liang Sun since Oct 2016

VarDiff compare two or more samples and identify the unique SNPs and small Indels

usage:
	VarDiff [options] -c <control1[,control2,...]> -m <mutant1[,mutant2,...]> -o <output file>
	
example:
	python VarDiff.py -c result/1AF_ngs-td7ddx4etd_S2_L001_R1_001.fastq.out   -m result/1AF_sample.fq1_sample.fq2.out -o result/test_unique.txt
    python VarDiff.py -c result/AF_ngs-0r7txc096a_S1_L001_R1_001.fastq.out   -m result/AF_ngs-td7ddx4etd_S2_L001_R1_001.fastq.out -o result/AF_unique.txt
	
"""

import sys
import os
import argparse
import time

#define mutant
class Mutant:
    def __init__(self, pos):
        self.name = pos
        self.chr = ""
        self.pos = int
        self.id = ""
        self.ref = ""
        self.alt = []
        self.af = float
    def averageAF(self,af):
        self.af = (self.af + af)/2
        
#read file
def readmFile(files):
    var_m = {} #chr_position -> freq
    for file in files:
        FileIN = open(file,"r")
        FileIN.readline()
        for line in FileIN:
            data = line.strip().split("\t")
            pos = (data[0], int(data[1]))
            if var_m.has_key(pos):
                var_m[pos].averageAF(float(data[5]))
                # need some work here to distinguish two allele at the same position
            else:
                var_m[pos] = Mutant(pos)
                var_m[pos].chr = data[0]
                var_m[pos].pos = int(data[1])
                var_m[pos].id = data[2]
                var_m[pos].ref = data[3]
                var_m[pos].alt = [data[4]]
                var_m[pos].af = float(data[5])
        FileIN.close()
    return var_m

def readcFile(files):
    var_c = {} #chr_position -> [allele]
    for file in files:
        FileIN = open(file,"r")
        FileIN.readline()
        for line in FileIN:
            data = line.strip().split("\t")
            pos = (data[0], int(data[1]))
            if var_c.has_key(pos):
                if data[4] in var_c[pos]:
                    continue
                else:
                    var_c[pos].append(data[4])
            else:
                var_c[pos] = [data[4]]
        
        FileIN.close()
    return var_c

#compare the mutant file and control file
def compareVar(var_m, var_c):
    var_u = {}
    for pos in var_m:
        if var_c.has_key(pos):
            continue
        else:
            var_u[pos] = var_m[pos]
            
    return var_u
        
#got the right format of input file from Yinbing

#output the right format which circos and accept

############################# declare variables #############################
t0 = time.time()

cfile = ""
mfile = ""
var_m = {} #chr_position -> freq
var_c = {} #chr_position -> allele
var_u = {} #chr_position -> freq   unique snp/indels

############################## process command line arguments #############################
parser = argparse.ArgumentParser(description="VarDiff compare two or more samples and identify the unique SNPs and Indels")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-c', action='store', dest='cfile', nargs='+', help="one or multiple wild type AF files")
parser.add_argument('-m', action='store', dest='mfile', nargs='+', help="one or multiple mutant sample AF files")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for unique AF files")
#parser.add_argument('-r', action='store', default='8', type=float, dest='cutoff', help="reads covrage cutoff")
args = parser.parse_args()

cfile = args.cfile
mfile = args.mfile
ofile = args.ofile
#cutoff = args.cutoff (how many )


if len(sys.argv) == 1:
	print parser.print_help()
	sys.exit("error: give me your input and output files!")

############################# read AF files #############################
else :
    var_m = readmFile(mfile)
    var_c = readcFile(cfile)

var_u = compareVar(var_m,var_c)

############################# write unique variance to output file #############################

ofile_snp = ofile+".snp"
ofile_del = ofile+".del"
FileOUT = open(ofile,"w")
FileOUT_snp = open(ofile_snp,"w")
FileOUT_del = open(ofile_del,"w")
FileOUT.write("CHROM\tPOS\tID\tREF\tALT\tAF_AVE\n")
FileOUT_snp.write("CHROM\tPOS\tID\tREF\tALT\tAF_AVE\n")
FileOUT_del.write("CHROM\tPOS_ST\tPOS_END\tDEL_LEN\tID\tREF\tALT\tAF_AVE\n")

for pos in sorted(var_u):
    FileOUT.write(var_u[pos].chr+"\t"+str(var_u[pos].pos)+"\t"+var_u[pos].id+"\t"+var_u[pos].ref+"\t"+str(var_u[pos].alt)+"\t"+ str(var_u[pos].af)+"\n")
    del_len = 0
    if len(var_u[pos].alt[0]) > len(var_u[pos].ref):
        #insertion
        continue
    elif len(var_u[pos].alt[0]) < len(var_u[pos].ref):
        #deletion
        del_len = len(var_u[pos].ref) - len(var_u[pos].alt[0])
        pos_end = int(var_u[pos].pos) + del_len
        FileOUT_del.write(var_u[pos].chr+"\t"+str(var_u[pos].pos)+"\t" + str(pos_end)+"\t"+str(del_len)+"\t"+var_u[pos].id+"\t"+var_u[pos].ref+"\t"+str(var_u[pos].alt)+"\t"+ str(var_u[pos].af)+"\n")
    else:
        FileOUT_snp.write(var_u[pos].chr+"\t"+str(var_u[pos].pos)+"\t"+var_u[pos].id+"\t"+var_u[pos].ref+"\t"+str(var_u[pos].alt)+"\t"+ str(var_u[pos].af)+"\n")
        
    
FileOUT.close()
FileOUT_snp.close()
FileOUT_del.close()
t1 = time.time()

totalTime = t1-t0
print "Run time: " + str(totalTime) + "sec"