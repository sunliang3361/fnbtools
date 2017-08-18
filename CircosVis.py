#!/usr/bin/env python

"""
to do:
1. check the number of big deletion, if too big, need to filter out smaller big deletions


"""

#encoding:utf-8
"""
CircosVis.py
Copyright (c) 2017 Noble Research Institute, LLC.

Created by Liang Sun on 2016-10-13
Updated and maintained by Liang Sun since Oct 2016

CircosVis visualizes small variance and large deletions

This portion of FNBTools depends on Circos, Copyright (c) 2004-2016 Martin Krzywinski, GPL License

usage:
	CircosVis [options] -c <control1[,control2,...]> -m <mutant1[,mutant2,...]> -o <output file>
	
example:
	python CircosVis.py  -l fnb/fnb.mt4_chr1_alldeletion_20x_annot.bed -o fnb_circos.png
"""

import sys
import os
import argparse
import time


        
def pickDeletionSize(count,lfile,minLen,maxLen):
    #minLen = 100
    while(count>25000):
        minLen = minLen + 10
        FileIN2 = open(lfile,"rU")
        FileIN2.readline()
        count = 0
        for line in FileIN2:
            data = line.strip().split("\t")
            if int(data[4])>=minLen and int(data[4])<=maxLen:
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
        if 'scaffold' in data[1]:
            continue
        if int(data[4]) < minLen:
            continue
        if int(data[4]) > minLen and int(data[4]) < maxLen:
            FileOUT.write(data[1]+"\t"+data[2]+"\t"+data[3]+"\t"+data[4]+"\n")

                
    FileIN.close()
    FileOUT.close()
def writeDelTextFile(trackfile,lfile,minLen,maxLen):
    print "*****Write text track to circos"
    FileOUT = open(trackfile,"w")
    FileIN = open(lfile,"rU")
    FileIN.readline()
    for line in FileIN:
        data = line.strip().split("\t")
        if 'scaffold' in data[1]:
            continue
        if minLen < int(data[4]) <= maxLen and data[11] != "[]":
            info = data[11].translate(None,"'[]").split(",")
            for id in info:
                FileOUT.write(data[1]+"\t"+data[2]+"\t"+data[3]+"\t"+id+"\n")
    FileIN.close()
    FileOUT.close()    
    
def processDel(lfile):
    FileIN1 = open(lfile,"rU")
    #count = sum(1 for line in FileIN1)
    FileIN1.readline()
    count1 = 0
    count100 = 0
    count1000 = 0
    for line in FileIN1:
        data = line.strip().split("\t")
        if int(data[3]) < 100:
            count1 = count1 + 1 
        elif int(data[3]) < 1000:
            count100 = count100 + 1
        else:
            count1000 = count1000 + 1
    FileIN1.close()
    # find the minimum length of deletion and wirte the deletion track files
    
    # minLen1 =  pickDeletionSize(count100,lfile,1,maxLen)
    # writeDelTrackFile("vis/DelTrack1.txt",lfile,minLen100,maxLen)

    minLen1 =  pickDeletionSize(count100,lfile,1,100)
    writeDelTrackFile("vis/DelTrack1.txt",lfile,minLen1,100)
    #write text file
    writeDelTextFile("vis/DelTrack1_text.txt",lfile,minLen1,100)
    
    minLen100 =  pickDeletionSize(count100,lfile,100,1000)
    writeDelTrackFile("vis/DelTrack100.txt",lfile,minLen100,1000)
    #write text file
    writeDelTextFile("vis/DelTrack100_text.txt",lfile,minLen100,1000)
    
    minLen1000 =  pickDeletionSize(count1000,lfile,1000,1500000)
    writeDelTrackFile("vis/DelTrack1000.txt",lfile,minLen1000,1500000)
    #write text file
    writeDelTextFile("vis/DelTrack1000_text.txt",lfile,minLen1000,1500000)
    

    


    
def createConf(ofile): #create in vis folder
    tempt_start = """
         # 1.2 IDEOGRAM LABELS, TICKS, AND MODULARIZING CONFIGURATION

        karyotype = vis/karyotype.mt.txt
        
        chromosomes_units = 1000000
        
        
        <plots>
        #text for dels from 1 to 100bp
        <plot>
        type       = text
        color      = 160, 0, 0
        label_font = condensed
        
        file = vis/DelTrack1_text.txt
        orientation = in
        
                r1   = 0.99r
                r0   = 0.80r
        
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
        #########histgram for dels from 1 to 100bp
        <plot>
        # The type sets the format of the track.
        
        type = histogram
        thickness = 2
        file = vis/DelTrack1.txt
        extend_bin = yes
        fill_color = 0, 0, 254
        color      = grey
        orientation = in
        r1   = 0.80r
        r0   = 0.65r
        
        # Histograms can have both a fill and outline. The default outline is 1px thick black. 
        <axes>
        <axis>
        spacing   = 0.4r
        color     = lgrey
        thickness = 1
        </axis>
        </axes>
        </plot>
        

        ######## text for dels from 100 to 1000bp
        <plot>
        type       = text
        color      = 0,150,0
        label_font = condensed
        
        file = vis/DelTrack100_text.txt
        orientation = in
        
                r1   = 0.65r
                r0   = 0.55r
        
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
        
        #########histgram for dels from 100 to 1000bp
        
        <plot>
        # The type sets the format of the track.
        
        type = histogram
        thickness = 2
        file = vis/DelTrack100.txt
        extend_bin = yes
        fill_color = 0, 200, 0
        color      = grey
        orientation = in
        r1   = 0.55r
        r0   = 0.40r
        # Histograms can have both a fill and outline. The default outline is 1px thick black. 
        <axes>
        <axis>
        spacing   = 0.4r
        color     = lgrey
        thickness = 1
        </axis>
        </axes>
        </plot>


        ######### text for dels from 1000 to 1,500,000bp
        <plot>
        type       = text
        color      = 0,0,190
        label_font = condensed
        
        file = vis/DelTrack1000_text.txt
        orientation = in
        
                r1   = 0.40r
                r0   = 0.30r
        
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
        
        #########histgram for dels from 1000 to 1,500,000bp
        <plot>
        # The type sets the format of the track.
        
        type = histogram
        thickness = 2
        file = vis/DelTrack1000.txt
        
        extend_bin = yes
        fill_color = 230,0,0
        color      = grey
        orientation = in
        r1   = 0.30r
        r0   = 0.15r
        # Histograms can have both a fill and outline. The default outline is 1px thick black. 
        <axes>
        <axis>
        spacing   = 0.4r
        color     = lgrey
        thickness = 1
        </axis>
        </axes>
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
#parser.add_argument('-s', action='store', dest='sfile', help="small variance file")
#parser.add_argument('-w', action='store', default= 50000, dest='wsize', help="windowsize for small variance file")
parser.add_argument('-l', action='store', dest='lfile', help="annoated large deletion file")
parser.add_argument('-o', action='store', dest='ofile', help="the output file name for circos image file")

args = parser.parse_args()

# sfile = args.sfile
# wsize = args.wsize
lfile = args.lfile
ofile = args.ofile

if len(sys.argv) == 1:
	print parser.print_help()
	sys.exit("error: give me your input and output files!")

############################# read AF files #############################
else :
    #processSmallVar(sfile,wsize)
    processDel(lfile)
    createConf(ofile)
    os.system("circos")



t1 = time.time()

totalTime = t1-t0
print "Total running time (sec): " + str(totalTime)