#!/usr/bin/env python


#encoding:utf-8
usage = """
VarDiff.py
created by Liang Sun on 2016-10-13
Copyright (c) 2016 Liang Sun. All rights reserved.
Updated and maintained by Liang Sun since Oct 2016

DelDiff compare two or more samples and identify the unique gap as the deletions

usage:
	DelDiff [options] -c <control1[,control2,...]> -m <mutant1[,mutant2,...]> -o <output file>
	
example:

	python DelDiff.py -c result/fnb.S1-FN3156_S1_R1.fastq.gz_gap.bedg result/fnb.S2-FN3164_S2_R1.fastq.gz_gap.bedg result/fnb.S4-FN6507_S4_R1.fastq.gz_gap.bedg -m result/fnb.S3-FN3196_S3_R1.fastq.gz_gap.bedg -f result/AF_unique_S3_del.txt -o result/fnb_deletion_S3.txt

"""

import sys
import os
import argparse
import time

sys.setrecursionlimit(100000000)


def read1File(delFile):
	print 'read gap file: '+ delFile
	FileIN = open(delFile, 'rU') or die ("can not open")
	chr_s={}
	for line in FileIN:
		data = line.strip().split("\t")
		length = int(data[2]) - int(data[1])
		# print data[2] + '\t' + data[1] + '\t***8'+ str(length)
		if length >= 1: #100
			if not chr_s.has_key(data[0]):
				chr_s[data[0]] = [(int(data[1]),int(data[2]))]
			else:
				chr_s[data[0]].append((int(data[1]),int(data[2])))
	FileIN.close()
	
	return chr_s


#read multiple inforamtive deletion files into one container
def readInfoDel(files,minCRR,minSMD,minInfo,minFR,type):
	chr_infodel = {}
	for f in files:
		#chr_infodel_c = read1File(f)
		#read informative reads, get 
		file = os.path.splitext(f)[0] + ".fix.bed"
		FileIN = open(file,"rU") or die ("can not open bed file!")
		#read gap file here also
		FileIN.readline()
		for line in FileIN:
			data = line.strip().split("\t")
			delInfo = "\t".join(data[4:])
			if data[2].find("[")>-1:
				continue
			clr_n = data[5].split(";")[0].split("=")[1]
			ccr_n = data[5].split(";")[1].split("=")[1]
			smd_n = data[5].split(";")[2].split("=")[1]
			flr_n = data[5].split(";")[3].split("=")[1]
			frr_n = data[5].split(";")[4].split("=")[1]
			if type == 'c':
				if chr_infodel.has_key(data[1]):
					#print data[1]+"---"+ data[2]+"------"+data[3]
					chr_infodel[data[1]].append((int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo))
					
				else:
					chr_infodel[data[1]] = [(int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo)]
			elif type == 'm':
				if int(clr_n) == 0 and int(smd_n) == 0 and int(ccr_n) < minCRR:          ####################cutoff to filter out inconfident deletions
					continue
				if int(clr_n) == 0 and int(smd_n) != 0 and int(smd_n) < minSMD:
					continue
				#if int(smd_n) != 0 and (int(clr_n)+int(smd_n)) < minInfo:
				if (int(clr_n)+int(smd_n)+int(ccr_n)) < minInfo:
					continue
				if int(flr_n) < minFR or int(frr_n) < minFR:
					continue
				if chr_infodel.has_key(data[1]):
					#print data[1]+"---"+ data[2]+"------"+data[3]
					chr_infodel[data[1]].append((int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo))
					
				else:
					chr_infodel[data[1]] = [(int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo)]			
	return chr_infodel 


def identifyHomo(chr_infodel_m,chr_gap_m,minDiff_rep,minDiff,orate):
	deletions = {}
	for chr in chr_infodel_m:
		for (brkpt,brkpt_end,delInfo) in chr_infodel_m[chr]:
			flag_gap_m = 0 #check the mutant gap files
			del_meta = ""
			gap_info = ""
			clr_n = int(delInfo.split("\t")[1].split(";")[0].split("=")[1])
			crr_n = int(delInfo.split("\t")[1].split(";")[1].split("=")[1])
			smd_n = int(delInfo.split("\t")[1].split(";")[2].split("=")[1])
			if not chr_gap_m or not chr_gap_m.has_key(chr):
				continue
			#check mutant gap file and decide whether it's homo
			if smd_n >= clr_n and smd_n > 0 : #small deletion
				for (gap_m_start,gap_m_end) in chr_gap_m[chr]:
					if abs(gap_m_start-brkpt) <= minDiff_rep: #3bp difference
						gap_info = str(gap_m_start)+"\t"+str(gap_m_end)
						flag_gap_m = 1
						break
			else: # big deletion
				flag_wiggle = 0
				gap_length_m = 0
				for (gap_m_start,gap_m_end) in chr_gap_m[chr]:
					#if (start >= brkpt and start<=brkpt_end) or (end>=brkpt and end <=brkpt_end): #or (start<=brkpt and end>=brkpt_end) not necessary
					if smd_n == 0 and clr_n == 0 : # only cross reads
						if gap_m_start-brkpt>=0 and gap_m_start-brkpt <= 50: #make the threshold loose
							#it is a homo deletion
							flag_wiggle = 1
							gap_info = str(gap_m_start)+"\t"+str(gap_m_end)
						if gap_m_start >= (brkpt-150) and gap_m_end <= (brkpt_end+150):
							gap_length_m = gap_length_m + (gap_m_end-gap_m_start)
					else:
						if gap_m_start-brkpt>=0 and gap_m_start-brkpt <= minDiff:
							#it is a homo deletion
							flag_wiggle = 1
							gap_info = str(gap_m_start)+"\t"+str(gap_m_end)
						if gap_m_start >= brkpt  and gap_m_end <= brkpt_end:
							gap_length_m = gap_length_m + (gap_m_end-gap_m_start)
						
				if flag_wiggle == 1 and float((gap_length_m)/float(brkpt_end-brkpt)) >= orate:
						flag_gap_m = 1
				#print chr+"-"+str(brkpt)+"-"+str(float((gap_length_m)/float(brkpt_end-brkpt)))
			if flag_gap_m == 1:
				#control informative deletion file has no deletion data or no deletion for this chromosome
				del_meta =  str(brkpt_end)+"\t"+delInfo+"\t"+gap_info + "\tYes\t-\t-"
				if deletions.has_key((chr,brkpt)):
					#print "Warniing! Two big deletions 2 at the same locations:"+ chr+"-"+str(brkpt)
					continue
				else:
					deletions[(chr,brkpt)] = del_meta
					
	return deletions

def filterDels(chr_infodel_m,chr_infodel_c,minDiff_rep,minDiff,type):
	deletions = {}
	for (chr,brkpt) in chr_infodel_m:
		flag_info_c = 0 #check the chr_infodel_c file
		brkpt_end = int(chr_infodel_m[(chr,brkpt)].split("\t")[0])
		
		# delInfo= chr_infodel_m[(chr,brkpt)].split("\t")[1]
		# gap_info = chr_infodel_m[(chr,brkpt)].split("\t")[2]
		clr_n = int(chr_infodel_m[(chr,brkpt)].split("\t")[2].split(";")[0].split("=")[1])
		crr_n = int(chr_infodel_m[(chr,brkpt)].split("\t")[2].split(";")[1].split("=")[1])
		smd_n = int(chr_infodel_m[(chr,brkpt)].split("\t")[2].split(";")[2].split("=")[1])
		if type == "all": # filter out all deletions
			if not chr_infodel_c or not chr_infodel_c.has_key(chr):
				pass
			else:
				if smd_n>= clr_n and smd_n > 0:
					for (brkpt_c,brkpt_end_c,delInfo_c) in chr_infodel_c[chr]:
						if abs(brkpt-brkpt_c) <= minDiff_rep:
							flag_info_c = 1
							break
				elif clr_n == 0 and smd_n == 0:
					for (brkpt_c,brkpt_end_c,delInfo_c) in chr_infodel_c[chr]:
						clr_n_c = int(delInfo_c.split("\t")[1].split(";")[0].split("=")[1])
						crr_n_c = int(delInfo_c.split("\t")[1].split(";")[1].split("=")[1])
						smd_n_c = int(delInfo_c.split("\t")[1].split(";")[2].split("=")[1])
						if clr_n_c == 0 and smd_n_c == 0:
							if abs(brkpt-brkpt_c) <= 150: #make the threshold loose
								flag_info_c = 1
								break
						else:
							if abs(brkpt-brkpt_c) <= minDiff:
								flag_info_c = 1
								break
				else:
					for (brkpt_c,brkpt_end_c,delInfo_c) in chr_infodel_c[chr]:
						if abs(brkpt-brkpt_c) <= minDiff:
							flag_info_c = 1
							break
		elif type == "homo": #filter out homo deletions from control samples
			#chr_infodel_c has different structure now and only contains homo deletions
			if not chr_infodel_c:
				pass
			if smd_n>= clr_n and smd_n > 0:
				#small deletions
				for (chr_c,brkpt_c) in chr_infodel_c:
					if abs(brkpt-brkpt_c) <= minDiff_rep:
						flag_info_c = 1
						break
			elif clr_n == 0 and smd_n == 0:
				for (chr_c,brkpt_c) in chr_infodel_c:
					clr_n_c = int(chr_infodel_c[(chr_c,brkpt_c)].split("\t")[2].split(";")[0].split("=")[1])
					crr_n_c = int(chr_infodel_c[(chr_c,brkpt_c)].split("\t")[2].split(";")[1].split("=")[1])
					smd_n_c = int(chr_infodel_c[(chr_c,brkpt_c)].split("\t")[2].split(";")[2].split("=")[1])
					if clr_n_c == 0 and smd_n_c == 0:
						if abs(brkpt-brkpt_c) <= 150: #make the threshold loose
							flag_info_c = 1
							break
					else:
						if abs(brkpt-brkpt_c) <= minDiff:
							flag_info_c = 1
							break
			else:
				for (chr_c,brkpt_c) in chr_infodel_c:
					if abs(brkpt-brkpt_c) <= minDiff:
						flag_info_c = 1
						break
		#deletions[(chr,brkpt)] =
		#print out final deletions
		#del_meta =  str(brkpt_end)+"\t"+delInfo+"\t"+gap_info+"\t"
		del_meta = "\t".join(chr_infodel_m[(chr,brkpt)].split("\t")[:5]) + "\t"
		if flag_info_c == 1:
			del_meta = del_meta + "Yes\tYes\tNo"
		else:
			del_meta  = del_meta + "Yes\tNo\tYes"

		if deletions.has_key((chr,brkpt)):
			#print "Warniing! Two big deletions 2 at the same locations:"+ chr+"-"+str(brkpt)
			continue
		else:
			deletions[(chr,brkpt)] = del_meta
			
	return deletions

############################# define all variables #############################
t0 = time.time()


parser = argparse.ArgumentParser(description="DelDiff compare two or more samples and identify the unique gap as the deletions")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-c', action='store', dest='cfile', nargs='*', help="one or multiple wild type bedg files, filter out homo and heter deletions")
parser.add_argument('-x', action='store', dest='xfile', nargs='*', help="one or multiple wild type bedg files, filter out homo deletions only")
parser.add_argument('-m', action='store', dest='mfile', nargs='*', help="one mutant sample bedg file")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for big file")
parser.add_argument('-r', action='store', default='0.9', type=float, dest='orate', help="gap overlapping rates")
parser.add_argument('-d', action='store', default='20', type=int, dest='minDiff', help="the minimal distance between the breakpoint and the start position of gap")
parser.add_argument('-b', action='store', default='3', type=int, dest='minCRR', help="the minimal crossed reads when there is no clipped reads [default:3]")
parser.add_argument('-s', action='store', default='3', type=int, dest='minSMD', help="the minimal small deletion reads [default:3]")
parser.add_argument('-i', action='store', default='3', type=int, dest='minInfo', help="the minimal total number of informative reads (clipped and small reads [default:2]")
parser.add_argument('-f', action='store', default='1', type=int, dest='minFR', help="the minimal flanking reads up and downstream of deletions")
parser.add_argument('-a', action='store', default='0', type=int, dest='allHomo', help="print all homo deletions in mutant including the deletions exist in control sample")
args = parser.parse_args()

cfile = ""
mfile = ""

chr_gap_c = {} # this is the all chromose and their deletion postion for control samples
chr_gap_m = {} # this is the all chromose and their deletion postion for mutant samples

chr_infodel_c = {} #chr_infodel[chr] = (breakpoint,support Reads)
chr_infodel_m = {}

deletions = {}
deletions_fix = {}

cfile = args.cfile
xfile = args.xfile
mfile = args.mfile
ofile = args.ofile
orate = args.orate #overlapping between control gap and the mutant inforamtive deletion
minDiff = args.minDiff #the minimal distance between the start position of informative deletion and the gap in the mutant.
minCRR = args.minCRR # if there is no clipped read, the minimal number of Crossed reads
minSMD = args.minSMD
minInfo = args.minInfo
minFR = args.minFR
allHomo = args.allHomo
minDiff_rep = 5 # the minimal distance difference to distinguish the informative deletion and small deletions

print "Use overlapping rate:"+ str(orate)
# print "cfile:"+str(len(cfile))
# print "mfile:"+str(len(mfile))
# print "orate:"+str(orate)	



if len(sys.argv) == 1:
	print parser.print_help()
	sys.exit("error: give me your input and output files!")

############################# more (>=1) samples vs more  (>=1) samples #############################
if cfile:
	chr_gap_m = read1File(mfile[0])
	
	file = mfile
	chr_infodel_c = readInfoDel(cfile,minCRR,minSMD,minInfo,minFR,'c') #here we have all homo and heter deletions
	chr_infodel_m = readInfoDel(file,minCRR,minSMD,minInfo,minFR,'m')
	chr_infodel_homo = identifyHomo(chr_infodel_m,chr_gap_m,minDiff_rep,minDiff,orate)
	deletions = filterDels(chr_infodel_homo,chr_infodel_c,minDiff_rep,minDiff,"all")
elif xfile:
	#print "only filter out homo deletions from control samples"
	chr_infodel_homo_m = {}
	chr_infodel_homo_c = {} #here we need all homo deletions

	file = mfile
	chr_gap_m_tmp = read1File(mfile[0])
	chr_infodel_m_tmp = readInfoDel(file,minCRR,minSMD,minInfo,minFR,'m')
	chr_infodel_homo_tmp = identifyHomo(chr_infodel_m_tmp,chr_gap_m_tmp,minDiff_rep,minDiff,orate)
	chr_infodel_homo_m.update(chr_infodel_homo_tmp)		
	#identify homo deletions in mutant samples:
	for f in xfile:
		file = [f]
		chr_gap_c_tmp = read1File(f)
		chr_infodel_c_tmp = readInfoDel(file,minCRR,minSMD,minInfo,minFR,'m')
		chr_infodel_homo_tmp = identifyHomo(chr_infodel_c_tmp,chr_gap_c_tmp,minDiff_rep,minDiff,orate)
		chr_infodel_homo_c.update(chr_infodel_homo_tmp)		
	#fitler out homo deletions in controls
	deletions = filterDels(chr_infodel_homo_m,chr_infodel_homo_c,minDiff_rep,minDiff,"homo")
else:
	print "no control file provided!"
	chr_gap_m = read1File(mfile[0])
	chr_infodel_m = readInfoDel(mfile,minCRR,minSMD,minInfo,minFR,'m')
	deletions = identifyHomo(chr_infodel_m,chr_gap_m,minDiff_rep,minDiff,orate)
	

				
			
FileOUT = open(ofile,'w') or die ("can not open file")

#FileFiltOUT = open(delFileFiltOUT,'w')
print ("Remove duplicate deletion!")

for (chr,brkpt) in deletions:
	#print "input:"+str((chr,brkpt))
	
	brkpt_end = int(deletions[(chr,brkpt)].split("\t")[0])
	clr_n = int(deletions[(chr,brkpt)].split("\t")[2].split(";")[0].split("=")[1])
	crr_n = int(deletions[(chr,brkpt)].split("\t")[2].split(";")[1].split("=")[1])
	smd_n = int(deletions[(chr,brkpt)].split("\t")[2].split(";")[2].split("=")[1])
	# print "input_end:"+str(brkpt_end)
	# print "clr_n:"+str(clr_n)
	# print "smd_n:" + str(smd_n)
	flag = 0
	for (chr1,brkpt1) in deletions:
		if chr1 == chr:
			brkpt_end1 = int(deletions[(chr1,brkpt1)].split("\t")[0])
			clr_n1 = int(deletions[(chr1,brkpt1)].split("\t")[2].split(";")[0].split("=")[1])
			crr_n1 = int(deletions[(chr1,brkpt1)].split("\t")[2].split(";")[1].split("=")[1])
			smd_n1 = int(deletions[(chr1,brkpt1)].split("\t")[2].split(";")[2].split("=")[1])
			if (brkpt >=brkpt1 and brkpt <= brkpt_end1) or (brkpt_end >=brkpt1 and brkpt_end <=brkpt_end1) or (brkpt<=brkpt1 and brkpt_end>=brkpt_end1):
				#print "why:"+str((chr1,brkpt1))+"--"+str(brkpt_end1)
				if (clr_n+smd_n+crr_n)<(clr_n1 + smd_n1+crr_n1):
					flag  = 1
				elif (clr_n+smd_n+crr_n) == (clr_n1 + smd_n1+crr_n1):
					if (clr_n+smd_n) < (clr_n1 + smd_n1):
						flag = 1
		else:
			continue
	if flag == 1:
		continue
	else:
		#print "output:"+str((chr,brkpt))
		deletions_fix[(chr,brkpt)] = deletions[(chr,brkpt)]
	

#FileOUT.write('Del#\tchr\tstart_position\tend_position\tdeletionLength\tbreakpoint_pos\tsupportRead\tdel_mutant\tdel_control\tHomo_Unique\n')
FileOUT.write('DEL#\tChr\tBreakpointStart\tBreakpointEnd\tDeletionLength\tSuppRead#\tGapStarts_position\tGapEnd_position\tDel_mutant\tDel_control\tHomo_Unique\n')

Del_num = 1;			
for (chr,brkpt) in sorted(deletions_fix):
	if allHomo == 0:
		if cfile or xfile:
			unique = deletions_fix[(chr,brkpt)].split("\t")[-1]
			if unique == "Yes":
				FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(brkpt)+"\t"+deletions_fix[(chr,brkpt)]+"\n")
				Del_num = Del_num + 1
		else:
			FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(brkpt)+"\t"+deletions_fix[(chr,brkpt)]+"\n")
			Del_num = Del_num + 1
	else:
		FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(brkpt)+"\t"+deletions_fix[(chr,brkpt)]+"\n")
		Del_num = Del_num + 1

			
			
				

FileOUT.close()

t1 = time.time()

totalTime = t1-t0
print "Total running time (sec): " + str(totalTime)