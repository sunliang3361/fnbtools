#!/usr/bin/env python


#encoding:utf-8
"""
VarAnnot.py
created by Liang Sun on 2016-10-13
Copyright (c) 2016 Liang Sun. All rights reserved.
Updated and maintained by Liang Sun since Oct 2016

VarAnnot annotate small variants (snp&small indel) and big deletions

usage:
	VarAnnot [options] -t type -i <input file> -o <output file>
	
example:
	python VarAnnot.py -i result/AF_unique.txt -t s -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/AF_unique_annot.txt
    #python VarAnnot.py  -i result/fnb_big_deletion.txt -t b -f Mtruncatula_285_Mt4.0v1.gene.gff3 -o result/fnb_big_deletion_annot.txt

"""

import sys
import os
import argparse
import time


        
#read file
def processGFF(file):
    geneInfo = {} # chr = [start, end, geneID]
    FileIN = open(file,"rU")
    for line in FileIN:
        if "##" in line:
            continue
        else:
            data = line.strip().split("\t")
            if data[2] == "gene":
                geneID = data[8].strip().split(";")[1].split("=")[1]
                if geneInfo.has_key(data[0]):
                    geneInfo[data[0]].append((int(data[3]),int(data[4]),geneID))
                else:
                    geneInfo[data[0]] = [(int(data[3]),int(data[4]),geneID)]
    
    FileIN.close()
    return geneInfo
        
def processVar(file):
    varInfo = {}
    FileIN = open(file,"rU")
    FileIN.readline()
    for line in FileIN:
        data = line.strip().split("\t")
        if data[1].find(":")>-1:
            id = data[1].split(":")[0]
            varInfo[(id,int(data[2]),int(data[3]))] = line.strip()
        else:
            varInfo[(data[1],int(data[2]),int(data[3]))] = line.strip()
    
    FileIN.close()
    return varInfo

def annotateVar(varInfo,geneInfo):  
    geneAnnot = {}
    for chr,start,end in varInfo:
        #print chr,start,end
        flag = 0
        if not geneInfo.has_key(chr):
            geneAnnot[(chr,start,end)] =[]
            continue
        for s,e,id  in geneInfo[chr]:
            if (start>= s and start <= e) or (end>=s and end <= e):
                if geneAnnot.has_key((chr,start,end)):
                    geneAnnot[(chr,start,end)].append(id)
                else:
                    flag = 1
                    geneAnnot[(chr,start,end)] = [id]
        if flag == 0:
            geneAnnot[(chr,start,end)] =[]
    
    return geneAnnot
                
            
            
        
    
############################# declare variables #############################
t0 = time.time()

ifile = ""
ffile = ""
ofile = ""
geneInfo = {} # chr = geneID ,start, end
varInfo = {} # (chr, pos) = info
geneAnnot = {} # (chr, pos) = geneID

############################## process command line arguments #############################
parser = argparse.ArgumentParser(description="VarAnnot annotate small variants (snp&small indel) and big deletions")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-i', action='store', dest='ifile', help="one unique small/big variance file")
parser.add_argument('-t', action='store', dest='ftype', help="b for big deletion")
parser.add_argument('-f', action='store', dest='ffile', help="gff3 annotation file")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for unique AF files")
#parser.add_argument('-r', action='store', default='8', type=float, dest='cutoff', help="reads covrage cutoff")
args = parser.parse_args()

ifile = args.ifile
ffile = args.ffile
ftype = args.ftype
ofile = args.ofile
#cutoff = args.cutoff (how many )


if len(sys.argv) == 1:
	print parser.print_help()
	sys.exit("error: give me your input and output files!")

############################# read AF files #############################
else :
    geneInfo = processGFF(ffile)
    varInfo = processVar(ifile)
    geneAnnot = annotateVar(varInfo,geneInfo)
    #print geneAnnot

############################# write unique variance to output file #############################
FileOUT = open(ofile,"w")
if ftype == "s":
    #FileOUT.write("#CHROM\tPOS\tID\tREF\tALT\tAF_AVE\tGENE\n")
    print "There is no s potion here"
elif ftype == "b":
    FileOUT.write("DEL#\tChr\tBreakpointStart\tBreakpointEnd\tDeletionLength\tSuppRead#\tGapStarts_position\tGapEnd_position\tDel_mutant\tDel_control\tHomo_Unique\tGenes\n")
for chr, start,end in sorted(varInfo):
    FileOUT.write(varInfo[(chr,start,end)]+"\t"+str(geneAnnot[(chr,start,end)])+"\n")
    
FileOUT.close()

t1 = time.time()
totalTime = t1-t0
print "Run time: " + str(totalTime) + "sec"
