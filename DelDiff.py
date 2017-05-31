#!/usr/bin/env python

# automatically detect the different deletion (zero coverage from bam file)
# if two deletions have less than 50% overlapping regions and both of them are >100bp, we will output it
# compare with the other three sample, if this deletion appears less than three time, we filter it out.

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

def findMin(g1,g2):
	#maxGap = []
	# print "g1:"
	# print g1
	# print "g2:"
	# print g2
	if not g1:
		return []
	elif not g2:
		return []
	else:
		gap1 = g1.pop(0)
		gap2 = g2.pop(0)
		s1 = gap1[0]
		e1 = gap1[1]
		s2 = gap2[0]
		e2 = gap2[1]
		if e1 <= s2:
			#maxGap.append((s1,e1))
			g2.insert(0,(s2,e2))
			#maxGap = maxGap + [(s1,e1)] + findMax(g1,g2)
			return findMin(g1,g2)
		elif e1 >= s2 and e1 <= e2:
			gap = (max(s1,s2),e1)	

			return [gap] + findMin(g1,g2)
		elif s1 <= s2 and e1 >= e2:
			return [(s2,e2)]+findMin(g1,g2)
		elif s1 >= s2 and s1 <= e2:
			gap = (s1,min(e1,e2))

			return [gap] + findMin(g1,g2)
		elif s1 >= e2:

			g1.insert(0,(s1,e1))
			#maxGap =maxGap + [(s2,e2)] + findMax(g1,g2)
			return findMin(g1,g2)


def overlapGap(chr_s1,chr_s2):
	chr_out = {}
	# print "chr_s1:"
	# print chr_s1
	# print "chr_s2:"
	# print chr_s2
	
	for chr in chr_s1:
		#print 'processing: '+chr
		if chr_s2.has_key(chr):
			chr_out[chr] = findMin(chr_s1[chr],chr_s2[chr])
		else:
			continue
			
	return chr_out

def findMax(g1,g2):
	#maxGap = []
	# print "g1:"
	# print g1
	# print "g2:"
	# print g2
	if not g1:
		#merge g2
		#maxGap = maxGap + g2
		#return maxGap

		return g2
	elif not g2:
		#merge g1
		#maxGap = maxGap + g1
		#return maxGap
		return g1
	else:
		gap1 = g1.pop(0)
		gap2 = g2.pop(0)
		s1 = gap1[0]
		e1 = gap1[1]
		s2 = gap2[0]
		e2 = gap2[1]
		if e1 <= s2:
			#maxGap.append((s1,e1))
			g2.insert(0,(s2,e2))
			#maxGap = maxGap + [(s1,e1)] + findMax(g1,g2)
			return [(s1,e1)] + findMax(g1,g2)
		elif e1 >= s2 and e1 <= e2:
			gap = (min(s1,s2),e2)	
			g2.insert(0,gap)
			#maxGap =maxGap + findMax(g1,g2)
			return findMax(g1,g2)
		elif s1 <= s2 and e1 >= e2:
			g1.insert(0,(s1,e1))
			return findMax(g1,g2)
		elif s1 >= s2 and s1 <= e2:
			gap = (s2,max(e1,e2))
			if e1>e2:
				g1.insert(0,gap)
			else:
				g2.insert(0,gap)
			#maxGap =maxGap + findMax(g1,g2)
			return findMax(g1,g2)
		elif s1 >= e2:
			#maxGap.append((s2,e2))
			g1.insert(0,(s1,e1))
			#maxGap =maxGap + [(s2,e2)] + findMax(g1,g2)
			return [(s2,e2)] + findMax(g1,g2)
			
def mergeGap(chr_s1,chr_s2):
	chr_out = {}

	for chr in chr_s1:
		#print 'processing: '+chr
		if chr_s2.has_key(chr):
			chr_out[chr] = findMax(chr_s1[chr],chr_s2[chr])
		else:
			chr_out[chr] = chr_s1[chr]
	for chr in chr_s2:
		if chr_s1.has_key(chr):
			continue
		else:
			chr_out[chr] = chr_s2[chr]
			
	return chr_out


def read_more_File(files, type):
	chr_s1 = {}
	flag = 0
	for f in files:
		f = f.strip()
		chr_s2 = {}
		if flag == 0:
			chr_s1 = read1File(f)
			flag = flag + 1
			continue
		else:
			chr_s2 = read1File(f)
		if type == "m":
			chr_out = {}
			chr_out = overlapGap(chr_s1,chr_s2)
			chr_s1.clear()
			chr_s1 = chr_out.copy()
		
		elif type == "c":
			chr_out = {}
			chr_out = mergeGap(chr_s1,chr_s2)
			chr_s1.clear()
			chr_s1 = chr_out.copy()
						
	return chr_s1

# def read_more_File(files):
# 	chr_gap = {}
# 	flag = 0
# 	for f in files:
# 		f = f.strip()
# 		print "Read gap file: "+f
# 		FileIN = open(f, 'rU') or die ("can not open")
# 		for line in FileIN:
# 			data = line.strip().split("\t")
# 			length = int(data[2]) - int(data[1])
# 			if length >= 1: #100
# 				if not chr_gap.has_key(data[0]):
# 					chr_gap[data[0]] = [(int(data[1]),int(data[2]))]
# 				else:
# 					chr_gap[data[0]].append((int(data[1]),int(data[2])))
# 		FileIN.close()
# 	return chr_gap

# def readSmVar(ffile):
# 	chr_filter = {}
# 	FileIN = open(ffile, 'rU') or die ("can not open filter file!")
# 	FileIN.readline()	
# 	for line in FileIN:
# 		data = line.strip().split("\t")
# 		if not chr_filter.has_key(data[0]):
# 			chr_filter[data[0]] = [(int(data[1]), int(data[2]))]
# 		else:
# 			chr_filter[data[0]].append((int(data[1]),int(data[2])))
# 	FileIN.close()
# 	return chr_filter


def readInfoDel(files,minCRR,minSMD,type):
	chr_infodel = {}
	for f in files:
		file = os.path.splitext(f)[0] + ".bed"
		FileIN = open(file,"rU") or die ("can not open bed file!")
		FileIN.readline()
		for line in FileIN:
			data = line.strip().split("\t")
			delInfo = "\t".join(data[4:])
			if data[2].find("[")>-1:
				continue
			clr_n = data[5].split(";")[0].split("=")[1]
			ccr_n = data[5].split(";")[1].split("=")[1]
			smd_n = data[5].split(";")[2].split("=")[1]
			if type == 'c':
				if chr_infodel.has_key(data[1]):
					#print data[1]+"---"+ data[2]+"------"+data[3]
					chr_infodel[data[1]].append((int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo))
					
				else:
					chr_infodel[data[1]] = [(int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo)]
			elif type == 'm':
				if int(clr_n) == 0 and int(smd_n) == 0 and int(ccr_n) < minCRR:          ####################cutoff to filter out inconfident deletions
					continue
				if int(smd_n) != 0 and (int(clr_n)+int(smd_n)) < minSMD:
					continue
				if chr_infodel.has_key(data[1]):
					#print data[1]+"---"+ data[2]+"------"+data[3]
					chr_infodel[data[1]].append((int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo))
					
				else:
					chr_infodel[data[1]] = [(int(float(data[2].replace(",",""))),int(float(data[3].replace(",",""))),delInfo)]
	return chr_infodel 
		



############################# define all variables #############################
t0 = time.time()

cfile = ""
mfile = ""

chr_gap_c = {} # this is the all chromose and their deletion postion for control samples
chr_gap_m = {} # this is the all chromose and their deletion postion for mutant samples

chr_smVar_c= {}
chr_smVar_m = {}
chr_infodel_c = {} #chr_infodel[chr] = (breakpoint,support Reads)
chr_infodel_m = {}


parser = argparse.ArgumentParser(description="DelDiff compare two or more samples and identify the unique gap as the deletions")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-c', action='store', dest='cfile', nargs='+', help="one or multiple wild type bedg files")
parser.add_argument('-m', action='store', dest='mfile', nargs='+', help="one or multiple mutant sample bedg files")
parser.add_argument('-f', action='store', dest='sfile', help="bcf del file used to filter small dels")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for big file")
parser.add_argument('-r', action='store', default='0.9', type=float, dest='orate', help="gap overlapping rates")
parser.add_argument('-d', action='store', default='20', type=int, dest='minDiff', help="the minimal distance between the breakpoint and the start position of gap")
parser.add_argument('-b', action='store', default='3', type=int, dest='minCRR', help="the minimal crossed reads when there is no clipped reads [default:3]")
parser.add_argument('-s', action='store', default='2', type=int, dest='minSMD', help="the minimal small deletion reads [default:2]")
parser.add_argument('-a', action='store', default='0', type=int, dest='allHomo', help="print all homo deletions in mutant including the deletions exist in control sample")
args = parser.parse_args()

cfile = args.cfile
mfile = args.mfile
sfile = args.sfile
ofile = args.ofile
orate = args.orate #overlapping between control gap and the mutant inforamtive deletion
minDiff = args.minDiff #the minimal distance between the start position of informative deletion and the gap in the mutant.
minCRR = args.minCRR # if there is no clipped read, the minimal number of Crossed reads
minSMD = args.minSMD
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
elif len(cfile) > 1 or len(mfile) > 1 :
	if len(cfile) > 1 :
		chr_gap_c = read_more_File(cfile,"c")  #find the merged gap regions
		#print "chr_gap_c:"
		#print chr_gap_c
	else:
		chr_gap_c = read1File(cfile[0])
		
	if len(mfile) > 1:
		chr_gap_m = read_more_File(mfile, "m")  #find the overlapping gap regions
		#print "chr_gap_m:"
		#print chr_gap_m
	else:
		chr_gap_m = read1File(mfile[0])


############################# one sample vs one sample #############################
elif len(cfile) == 1 and len(mfile) == 1:
	chr_gap_c = read1File(cfile[0])

	chr_gap_m = read1File(mfile[0])
	
# elif len(cfile) >= 1 and len(mfile) >= 1 :
# 	chr_gap_c = read_more_File(cfile)  #find the merged gap regions
# 	chr_gap_m = read_more_File(mfile)  #find the overlapping gap regions
	
chr_infodel_c = readInfoDel(cfile,minCRR,minSMD,'c')
chr_infodel_m = readInfoDel(mfile,minCRR,minSMD,'m')

# sfile_ctrl_del = sfile+".ctrl.del"
# sfile_mut_del = sfile+".mut.del"
# chr_smVar_c = readSmVar(sfile_ctrl_del)
# chr_smVar_m= readSmVar(sfile_mut_del)
	

### calculate the reliability of each identified 
#calculate by the overlapping between normal and mutant samples

#delFileOUT = 'result/significant_fnb.txt'

#need to figure out the cutoff??????????????????????????????
#delFileFiltOUT = 'result/significant_filt_fnb.txt'
#print chr_gap_c['chr2']

deletions = {}
deletions_fix = {}

for chr in chr_infodel_m:
	for (brkpt,brkpt_end,delInfo) in chr_infodel_m[chr]:
		flag_gap_m = 0 #check the mutant gap files
		flag_info_c = 0 #check the chr_infodel_c file
		flag_gap_c = 0 #check the control gap file
		
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
			#deletions[(chr,brkpt)] =
			del_meta =  str(brkpt_end)+"\t"+delInfo+"\t"+gap_info+"\t"
			if flag_info_c == 1:
				del_meta = del_meta + "Yes\tYes\tNo"
			else:
				del_meta  = del_meta + "Yes\tNo\tYes"
				#check whether the control gap file has gap at these area
				# gap_length_c = 0
				# if not chr_gap_c or not chr_gap_c.has_key(chr):
				# 	del_meta  = del_meta + "Yes\tNo\tYes"
				# else:
					# if smd_n: #small deletion
					# 	for (gap_c_start,gap_c_end) in chr_gap_c[chr]:
					# 		if gap_c_start <= brkpt and gap_c_end >= brkpt_end:
					# 			flag_gap_c = 1
					# 			break
					# 		elif (gap_c_start >= brkpt and gap_c_start <=brkpt_end) or (gap_c_end >=brkpt and gap_c_end <= brkpt_end):
					# 			#if float((min(gap_c_end,brkpt_end)-max(gap_c_start,brkpt))/float(brkpt_end-brkpt)) >= 0.5:
					# 			flag_gap_c = 1
					# 			break
					# 	if flag_gap_c == 1:
					# 		del_meta = del_meta + "Yes\tMaybe\tMaybe"
					# 	else:
					# 		del_meta = del_meta + "Yes\tNo\tYes"
					# else: # large deletion
					# 	#for (start_c,end_c) in chr_gap_c[chr]:
					# 	for (gap_c_start,gap_c_end) in chr_gap_c[chr]:
					# 		if gap_c_start <= brkpt and gap_c_end >= brkpt_end:
					# 			flag_gap_c = 1
					# 			break
					# 		elif (gap_c_start >= brkpt and gap_c_start <=brkpt_end) or (gap_c_end >=brkpt and gap_c_end<=brkpt_end):
					# 			#if float((min(end_c,brkpt_end)-max(start_c,brkpt))/float(brkpt_end-brkpt)) >= (1 - orate):
					# 			gap_length_c = gap_length_c + (min(gap_c_end,brkpt_end)-max(gap_c_start,brkpt))
					# 	if flag_gap_c == 1 or float(gap_length_c/float(brkpt_end-brkpt)) >= orate:
					# 		del_meta = del_meta + "Yes\tMaybe\tMaybe"
					# 	elif float(gap_length_c/float(brkpt_end-brkpt)) < float(1- orate): # 90% gap is not overlap with deletion in mutant 
					# 		del_meta = del_meta + "Yes\tNo\tYes"
					# 	else:
					# 		continue
					

			if deletions.has_key((chr,brkpt)):
				#print "Warniing! Two big deletions 2 at the same locations:"+ chr+"-"+str(brkpt)
				continue
			else:
				deletions[(chr,brkpt)] = del_meta

							
					
			
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
#DEL#\tChr\tBreakpointStart\tBreakpointEnd\tDeletionLength\tSuppRead#
#FileFiltOUT.write('chr\tstart_position\tend_position\tdeletionLength\tfreq\n')
Del_num = 1;			
for (chr,brkpt) in sorted(deletions_fix):
	if allHomo == 0:
		unique = deletions_fix[(chr,brkpt)].split("\t")[-1]
		if unique == "Yes":
			FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(brkpt)+"\t"+deletions_fix[(chr,brkpt)]+"\n")
			Del_num = Del_num + 1
	else:
		FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(brkpt)+"\t"+deletions_fix[(chr,brkpt)]+"\n")
		Del_num = Del_num + 1

			
			
#DEL#\tChr\tBreakpointStart\tBreakpointEnd\tDeletionLength\tSuppRead#\tGapStartstart_position\tGapEndend_position\tDel_mutant\tDel_control\tHomo_Unique
				

FileOUT.close()

t1 = time.time()

totalTime = t1-t0
print "Total running time (sec): " + str(totalTime)