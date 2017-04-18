#!/usr/bin/env python

# automatically detect the different deletion (zero coverage from bam file)
# if two deletions have less than 50% overlapping regions and both of them are >100bp, we will output it
# compare with the other three sample, if this deletion appears less than three time, we filter it out.

#encoding:utf-8
"""
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

	
	
def compareDel(chr_s1,chr_s2,orate):
	#chr_s1 is the mutanted sample
	chr_sig = {}
	for chr in chr_s1:
		# print chr_s1[chr]
		for s1, e1 in chr_s1[chr]:
			# print str(s1)+'------' + str(e1)
			if not chr in chr_s2:
				continue
			
			flag = 0
			for s2, e2 in chr_s2[chr]:
				if s2>=s1 and s2<=e1:
					flag = 1
					if e2< e1:
						olp = float(float(e2-s2)/float(e1-s1))
						#print str(s1)+'--'+str(e1)+'--'+str(s2)+'--'+str(e2)
						if olp < orate:
							# print chr+"\t"+str(s1)+"\t"+str(e1) +'--------' +chr+"\t"+str(s2)+"\t"+str(e2) 
							chr_sig[(chr,s1,e1)] = "low"
							break

					else:
						olp = float(float(e1-s2)/float(e1-s1))
						if olp < orate:
							# print chr+"\t"+str(s1)+"\t"+str(e1) +'--------' +chr+"\t"+str(s2)+"\t"+str(e2)
							chr_sig[(chr,s1,e1)] = "low"
							break
							# print str(s1)+'**'+str(e1)+'**'+str(s2)+'**'+str(e2)
							# print str(olp)
							# print "&&&&&&&&"
							#FileOUT.write(chr + '\t'+ str(s1)+ '\t'+ str(e1)+ '\t'+ str((e1-s1))+'\t'+ str(s2)+ '\t'+ str(e2)+ '\t'+str((e2-s2))+"\n")
				elif s2<=s1 and e2>=s1:
					flag = 1
					if e2<e1:
						olp = float(float(e2-s1)/float(e1-s1))
						if olp < orate:
							# print chr+"\t"+str(s1)+"\t"+str(e1) +'--------' +chr+"\t"+str(s2)+"\t"+str(e2)
							chr_sig[(chr,s1,e1)] = "low"
							break
							#FileOUT.write(chr + '\t'+ str(s1)+ '\t'+ str(e1)+ '\t'+ str((e1-s1))+'\t'+ str(s2)+ '\t'+ str(e2)+ '\t'+str((e2-s2))+"\n")
					else:
						break
							#FileOUT.write(chr + '\t'+ str(s1)+ '\t'+ str(e1)+ '\t'+ str((e1-s1))+'\t'+ str(s2)+ '\t'+ str(e2)+ '\t'+str((e2-s2))+"\n")

			if flag == 0:
				# print chr+"\t"+str(s1)+"\t"+str(e1) +'--------' +chr+"\t"+str(s2)+"\t"+str(e2)
				chr_sig[(chr,s1,e1)] = "high"					
				#FileOUT.write(chr + '\t'+ str(s1)+ '\t'+ str(e1)+ '\t'+ str((e1-s1))+'\tNA\tNA\tNA\n')
	return chr_sig
def readfilter(ffilter):
	chr_filter = {}
	FileIN = open(ffile, 'rU') or die ("can not open filter file!")
	FileIN.readline()
	global del_bcf_max
	
	for line in FileIN:
		data = line.strip().split("\t")
		del_bcf_max = max(del_bcf_max,int(data[3]))
		if not chr_filter.has_key(data[0]):
			chr_filter[data[0]] = [(int(data[1]), int(data[2]))]
		else:
			chr_filter[data[0]].append((int(data[1]),int(data[2])))
	FileIN.close()
	return chr_filter

def readInfoDel(files):
	chr_infodel = {}
	for f in files:
		file = os.path.splitext(f)[0] + ".bed"
		FileIN = open(file,"rU") or die ("can not open bed file!")
		FileIN.readline()
		for line in FileIN:
			data = line.strip().split("\t")
			if data[2].find("[")>-1:
				continue
			if chr_infodel.has_key(data[1]):
				#print data[2]+"------"+data[3]
				chr_infodel[data[1]].append((int(float(data[2].replace(",",""))),data[3]))
			else:
				chr_infodel[data[1]] = [(int(float(data[2].replace(",",""))),data[3])]
	return chr_infodel 
		



############################# define all variables #############################
t0 = time.time()

cfile = ""
mfile = ""
del_bcf_max = 0
chr_gap_c = {} # this is the all chromose and their deletion postion for control samples
chr_gap_m = {} # this is the all chromose and their deletion postion for mutant samples

chr_sig = {}
chr_filter = {}
chr_infodel_c = {} #chr_infodel[chr] = (breakpoint,support Reads)
chr_infodel_m = {}


parser = argparse.ArgumentParser(description="DelDiff compare two or more samples and identify the unique gap as the deletions")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-c', action='store', dest='cfile', nargs='+', help="one or multiple wild type bedg files")
parser.add_argument('-m', action='store', dest='mfile', nargs='+', help="one or multiple mutant sample bedg files")
parser.add_argument('-f', action='store', dest='ffile', help="bcf del file used to filter small dels")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for big file")
parser.add_argument('-r', action='store', default='0.5', type=float, dest='orate', help="gap overlapping rates")
args = parser.parse_args()

cfile = args.cfile
mfile = args.mfile
ffile = args.ffile
ofile = args.ofile
orate = args.orate

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


##################### compare the gaps between merged gap #####################
chr_sig = compareDel(chr_gap_m,chr_gap_c,orate)
chr_filter = readfilter(ffile)

chr_infodel_c = readInfoDel(cfile)
chr_infodel_m = readInfoDel(mfile)

### calculate the reliability of each identified 
#calculate by the overlapping between normal and mutant samples

#delFileOUT = 'result/significant_fnb.txt'

#need to figure out the cutoff??????????????????????????????
#delFileFiltOUT = 'result/significant_filt_fnb.txt'
FileOUT = open(ofile,'w') or die ("can not open file")

#FileFiltOUT = open(delFileFiltOUT,'w')

FileOUT.write('Del#\tchr\tstart_position\tend_position\tdeletionLength\tbreakpoint_pos\tsupportRead\tdel_mutant\tdel_control\tHomo_Unique\n')
#FileFiltOUT.write('chr\tstart_position\tend_position\tdeletionLength\tfreq\n')
Del_num = 1;
for (chr,start,end) in sorted(chr_sig):
	delsize = end-start
	#print "the biggest dels in bcf"+str(del_bcf_max)
	if delsize > del_bcf_max:
		flag = 0
		flag_c = 0
		if not chr_infodel_m:
			FileOUT.write(str(Del_num)+chr+"\t"+str(start)+"\t"+str(end)+"\t"+str((end-start))+"\t-\t-\t-\t-\t-\n")
			
		if not chr_infodel_m.has_key(chr):
			continue
		for (brkpt,sprd) in chr_infodel_m[chr]:
			if abs(start-brkpt)<=100:
				flag = 1
				break
		if not chr_infodel_c or not chr_infodel_c.has_key(chr):
			if flag ==1:
				FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(start)+"\t"+str(end)+"\t"+str((end-start))+"\t"+str(brkpt)+"\t"+sprd+"\tYes\tNo\tYes\n")
				Del_num = Del_num + 1
			continue
		for (brkpt_c,sprd_c) in chr_infodel_c[chr]:
			if abs(start-brkpt_c)<=50:
				flag_c = 1
		if flag == 1:
			FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(start)+"\t"+str(end)+"\t"+str((end-start))+"\t"+str(brkpt)+"\t"+sprd+"\tYes\t")
			Del_num = Del_num + 1
			if flag_c == 1:
				FileOUT.write("Yes\tNo\n")
			else:
				FileOUT.write("No\tYes\n")
		#FileOUT.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str((end-start))+"\t"+str(chr_sig[(chr,start,end)])+"\n")
	else:
		flag = 0
		if not chr_filter.has_key(chr):
			continue
		for (st, ed) in chr_filter[chr]:
			if (start >=st and start <=ed) or (end >=st and end <=ed):
				flag = 1
				break
		if flag == 1:
			FileOUT.write(str(Del_num)+"\t"+chr+"\t"+str(start)+"\t"+str(end)+"\t"+str((end-start))+"\t-\t-\tYes\tNO\tYes\n")
			Del_num = Del_num + 1
		else:
			continue
				
				
			

FileOUT.close()

t1 = time.time()

totalTime = t1-t0
print "Total running time (sec): " + str(totalTime)