#!/usr/bin/env python

"""
to do:
1. check the number of big deletion, if too big, need to filter out smaller big deletions


"""

#encoding:utf-8
"""
CircosVis.py
created by Liang Sun on 2016-10-13
Copyright (c) 2016 Liang Sun. All rights reserved.
Updated and maintained by Liang Sun since Oct 2016

CircosVis visualize small variance and large deletions

usage:
	VarDiff [options] -c <control1[,control2,...]> -m <mutant1[,mutant2,...]> -o <output file>
	
example:
	python CircosVis.py -s result/AF_unique.txt  -l result/fnb_big_deletion_annot.txt -o myfile.png
	python CircosVis.py -s result/AF_unique_S2.txt  -l result/fnb_big_deletion_S2_annot.txt -o myfile.png
"""

import sys
import os
import argparse
import time

# class VarInfo:
#     def __init__(self):
#         self.chr = ""
#         self.start = int
#         self.af = []
#         self.af_ave = float
#         
#     def averageAF():
#         return sum(self.af)/len(self.af)

#read file
def processSmallVar(sfile,wsize):
    chr_win = {} #chr_win[(chr,st,st+wsize)] = [AF]
    #del_hom = {} #del_hom[(chr,st,st+wsize)] = AF
    chr = []
    st = 0
    FileIN = open(sfile,"rU")
    FileIN.readline()
    for line in FileIN:
        data = line.strip().split("\t")
        if 'scaffold' in data[0]:
            continue
        if data[0] in chr:
            if int(data[1]) <= (st + wsize):
                
                if chr_win.has_key((data[0],st,(st + wsize))):
                    chr_win[(data[0],st,(st + wsize))].append(float(data[5]))
                else: #not nessariry
                    chr_win[(data[0],st,(st + wsize))] = [float(data[5])]
            else:
                st = st + wsize
                chr_win[(data[0],st,(st+wsize))] = [float(data[5])]
        else:
            chr.append(data[0])
            st = 0
            chr_win[(data[0],st,(st+wsize))] = [float(data[5])]
                
    FileIN.close()
    # check whether too much data to show in circos
    if len(chr_win) > 25000:
        exit("Error message: too many data to show in circos, you need to increase the window size parameter!")

    FileOUT = open('vis/VarTrack.txt','w')
    for chr,start,end in sorted(chr_win):
        af_ave = sum(chr_win[(chr,start,end)])/len(chr_win[(chr,start,end)])
        FileOUT.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(af_ave)+"\n")
    FileOUT.close()
    
    # FileOUT1 = open('vis/Del0Track.txt','w')
    # for chr,start,end in sorted(chr_win):
    #     af_ave = sum(chr_win[(chr,start,end)])/len(chr_win[(chr,start,end)])
    #     FileOUT.write(chr+"\t"+str(start)+"\t"+str(end)+"\t"+str(af_ave)+"\n")
    # FileOUT1.close()
        
def pickDeletionSize(count,lfile,minLen,maxLen):
    #minLen = 100
    while(count>25000):
        minLen = minLen + 10
        FileIN2 = open(lfile,"rU")
        FileIN2.readline()
        count = 0
        for line in FileIN2:
            data = line.strip().split("\t")
            if int(data[3])>=minLen and int(data[3])<=maxLen:
                count = count + 1
        FileIN2.close()
    return minLen

def writeDelTrackFile(trackfile,lfile,minLen,maxLen):
    print "*****The minimum Length of big deletions to be visualized in track" + trackfile + ": "+str(minLen) + "bp"
    FileOUT = open(trackfile,"w")
    FileIN = open(lfile,"rU")
    FileIN.readline()
    for line in FileIN:
        data = line.strip().split("\t")
        if 'scaffold' in data[0]:
            continue
        if int(data[3]) < minLen:
            continue
        if int(data[3]) > minLen and int(data[3]) < maxLen:
            if data[4] == "low":
                FileOUT.write(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+data[3]+"\t"+"fill_color=blue"+"\n")
            elif data[4] == "high":
                FileOUT.write(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+data[3]+"\t"+"fill_color=red" +"\n")
                
    FileIN.close()
    FileOUT.close()
def writeDelTextFile(trackfile,lfile,maxLen):
    print "*****Write text track to circos"
    FileOUT = open(trackfile,"w")
    FileIN = open(lfile,"rU")
    FileIN.readline()
    for line in FileIN:
        data = line.strip().split("\t")
        if 'scaffold' in data[0]:
            continue
        if int(data[3]) < maxLen or data[5] == "[]":
            continue
        else:
            info = data[5].translate(None,"'[]").split(",")
            for id in info:
                FileOUT.write(data[0]+"\t"+data[1]+"\t"+data[2]+"\t"+id+"\n")
    FileIN.close()
    FileOUT.close()    
    
def processLargeVar(lfile,maxLen):
    FileIN1 = open(lfile,"rU")
    #count = sum(1 for line in FileIN1)
    FileIN1.readline()
    count100 = 0
    count1000 = 0
    for line in FileIN1:
        data = line.strip().split("\t")
        if int(data[3]) < maxLen:
            count100 = count100 + 1
        else:
            count1000 = count1000 + 1
    FileIN1.close()
    # find the minimum length of deletion and wirte the deletion track files
    minLen100 =  pickDeletionSize(count100,lfile,100,maxLen)
    writeDelTrackFile("vis/DelTrack100.txt",lfile,minLen100,maxLen)
    
    minLen1000 =  pickDeletionSize(count1000,lfile,maxLen,10000000)
    writeDelTrackFile("vis/DelTrack1000.txt",lfile,minLen1000,10000000)
    
    #write text file
    writeDelTextFile("vis/DelTrack1000_text.txt",lfile,maxLen)
    

    


    
def createConf(ofile): #create in vis folder
    tempt_start = """
         # 1.2 IDEOGRAM LABELS, TICKS, AND MODULARIZING CONFIGURATION

        karyotype = vis/karyotype.mt.txt
        
        chromosomes_units = 1000000
        
        
        <plots>

        
        <plot>
        # The type sets the format of the track.
        
        type = histogram
        thickness = 2
        
        file = vis/VarTrack.txt
        extend_bin = yes
        fill_color = vdgrey
        color      = black
        
        r1   = 0.95r
        r0   = 0.80r
        
        <rules>
        
        <rule>
        condition = var(value) == 1.0
        color = red
        </rule>

        </rules>
        
        # Histograms can have both a fill and outline. The default outline is 1px thick black. 
        
        
        <axes>
        <axis>
        spacing   = 0.1r
        color     = lgrey
        thickness = 2
        </axis>
        </axes>
        </plot>
        
        <plot>
        # The type sets the format of the track.
        
        type = histogram
        thickness = 2
        file = vis/DelTrack100.txt
        #extend_bin = yes
        fill_color = vdgrey
        color      = grey
        
        r1   = 0.75r
        r0   = 0.60r
        

        
        # Histograms can have both a fill and outline. The default outline is 1px thick black. 
        <axes>
        <axis>
        spacing   = 0.1r
        color     = lgrey
        thickness = 2
        </axis>
        </axes>
        </plot>
        
        <plot>
        # The type sets the format of the track.
        
        type = histogram
        thickness = 2
        file = vis/DelTrack1000.txt
        
        #extend_bin = yes
        fill_color = vdgrey
        color      = grey
        
        r1   = 0.55r
        r0   = 0.40r
        

        
        # Histograms can have both a fill and outline. The default outline is 1px thick black. 
        <axes>
        <axis>
        spacing   = 0.1r
        color     = lgrey
        thickness = 2
        </axis>
        </axes>
        </plot>


<plot>
type       = text
color      = black
label_font = condensed

file = vis/DelTrack1000_text.txt
orientation = in

        r1   = 0.39r
        r0   = 0.20r

label_size = 12p

show_links     = yes
link_dims      = 2p,2p,4p,2p,2p
link_thickness = 2p
link_color     = red

label_snuggle         = yes
max_snuggle_distance  = 1r
snuggle_tolerance     = 0.25r
snuggle_sampling      = 2


</plot>
</plots>


        
        # don't include links
        
        <<include etc/ideogram.conf>>
        
        <<include etc/ticks.conf>>
        <image>
        <<include etc/image.conf>>
    """
    fileName = "file* = vis/"+ofile

    tempt_end = """                  
        </image>
        
        <<include etc/colors_fonts_patterns.conf>> 
        
        <<include etc/housekeeping.conf>> 
        
        data_out_of_range* = trim
    
    """
    FileOUT = open("circos.conf","w")
    FileOUT.write(tempt_start+"\n")
    FileOUT.write(fileName+"\n")
    FileOUT.write(tempt_end+"\n")
    FileOUT.close()
    
#got the right format of input file from Yinbing

#output the right format which circos and accept

############################# declare variables #############################
t0 = time.time()
sfile = ""
wfile = int
lfile = ""
ofile = ""

############################# check the availability of pre-required software #############################
#check circos module installed or not????????????????????????????????????????????????????
#

############################## process command line arguments #############################
parser = argparse.ArgumentParser(description="VarAnnot annotate small variants (snp&small indel) and big deletions")
parser.add_argument('--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-s', action='store', dest='sfile', help="small variance file")
parser.add_argument('-w', action='store', default= 50000, dest='wsize', help="windowsize for small variance file")
parser.add_argument('-l', action='store', dest='lfile', help="annoated large deletion file")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for circos image file")

args = parser.parse_args()

sfile = args.sfile
wsize = args.wsize
lfile = args.lfile
ofile = args.ofile

if len(sys.argv) == 1:
	print parser.print_help()
	sys.exit("error: give me your input and output files!")

############################# read AF files #############################
else :
    processSmallVar(sfile,wsize)
    processLargeVar(lfile,1000)
    createConf(ofile)
    os.system("circos")



t1 = time.time()

totalTime = t1-t0
print "Total running time (sec): " + str(totalTime)
