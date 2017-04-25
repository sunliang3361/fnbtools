# Installation

## Install from Source
1. git clone --- (or download fnbtools.tar.gz from Releases and run `tar -xzvf fnbtools.tar.gz`)
2. sudo ./Install.sh

## Install in Docker container
1. git clone ---
2. 


# Using FNBTools

### Step 1 - Align all reads to reference genome
Usage:
	fnbalign -n fnb -g genome.fa -1 control_forward.fq mutant_forward.fq -2 control_reverse.fq mutant_reverse.fq

Parameters:
	REQUIRED -g the genome sequence file in fasta format  
	REQUIRED -1 the paired read file 1
	REQUIRED -2 the paired read file 2
	REQUIRED -n the the name of you project

Example:
perl fnbalign.pl -n fnb -g mt4_chr1_2Mb.fa -1 mt4_chr1_raw_20x1.fq mt4_chr1_mut_20x1.fq -2 mt4_chr1_raw_20x2.fq mt4_chr1_mut_20x2.fq

### Step2 - Identify deletions in mutant sample(s)
Usage:
	fnbscan -n fnb -c control.bed -m mutant.bed
Parameters:
	REQUIRED -n the project name
	REQUIRED -c the control bed file
	REQUIRED -m the mutant bed file
	REQUIRED -o output file name for identified deletions

	OPTIONS:
	-f input annotation file to annotate deletions     
	-h print this help message 
	
Example:
perl fnbscan.pl -n fnb -c fnb/fnb.mt4_chr1_raw_20x1.bedg -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -f Mtruncatula_285_Mt4.0v1.gene.gff3

### Step3 - Visualize all identified deletions via Circos
Usage:
	CircosVis [options] -s <snp.vcf> -l <bigDeletion.vcf> -o <output file>
Parameters:
	REQUIRED -s the unique SNP file  
	REQUIRED -l the large deletion file
	REQUIRED -o the output image file
Example:
	python CircosVis.py -s result/AF_unique.txt  -l result/fnb_big_deletion_annot.txt -o myfile.png
