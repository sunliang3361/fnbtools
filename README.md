# Installation

## Install from Source
1. Run: `git clone ---` (or download fnbtools.tar.gz from Releases and run `tar -xzvf fnbtools.tar.gz`)
2. Run: `sudo ./Install.sh`
3. If prompted with "Would you like to configure as much as possible automatically? [yes]", type **yes** and then press Enter. This will automatically configure Perl's CPAN utility so that additional Perl modules can be installed.
4. Run: `. ~/.bashrc`

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

Output:
1. BAM file
There are one bam file for each sample in the Temp folder 
e.g. fnb.mt4_chr1_mut_20x1.ref.sort.bam, fnb.mt4_chr1_raw_20x1.ref.sort.bam
2. BEDG file
One Gap file for each sample in the project folder (e.g. fnb folder)
e.g. fnb.mt4_chr1_mut_20x1.bedg, fnb.mt4_chr1_raw_20x1.bedg
2. BED file
All preliminary deletion data is stored here 
There is one bed file for mutant samples and one bed file for control samples.
e.g mutant samples: fnb.mt4_chr1_mut_20x1.bed
	control samples: fnb.mt4_chr1_raw_20x1.bed
### Step 2 - Identify deletions in mutant sample(s)
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
fnbscan -n fnb -c fnb/fnb.mt4_chr1_raw_20x1.bedg -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -f Mtruncatula_285_Mt4.0v1.gene.gff3

Output:
1. BED file
All unique deletions will be output in this file under the project folder (User needs to give the project name by using Parameter -n, e.g. fnb).
e.g. fnb.mt4_chr1_alldeletion_20x.bed
DEL#: deletion number. All deletions will be sorted by the chromosome number and start postions
Chr: Chromosome number.	
BreakpointStart: Start position of deletions
BreakpointEnd: End position of deletions
DeletionLength: Deletion length
SuppRead#: CLR represent the clipped reads number; CRR represents discordant reads number; SMD represents small deletion reads number 
GapStarts_position: Start position of gaps	
GapEnd_position	: End position of gaps
Del_mutant: 'Yes' represents deletions found in mutant samples	
Del_control: 'Yes' represents deletions found in control samples; 'No' represents no deletions found in control samples
Homo_Unique: 'Yes' represents deletions only exist in mutant samples but not in control sample; 'No' represents deletions exist in both mutant and control sample
Genes (optional): This column will only show genes if deletions cover

2. BED file with annotation (option)
### Step3 - Visualize all identified deletions via Circos
Usage:
	CircosVis [options] -l <bigDeletion.vcf> -o <output file>
Parameters:
	REQUIRED -l the large deletion file
	REQUIRED -o the output image file
Example:
	circosvis  -l example/wen.S1_alldeletion_uniq_annot.bed -o fnb_circos.png
	circosvis  -l fnb/fnb.mt4_chr1_alldeletion_20x_annot.bed -o fnb_circos.png

Output:
The output circos image will be in vis folder with the name parameter 'o' you provided
If no error message reported, you can find the output circos image in vis folder. Since our sample data (fnb.mt4_chr1_alldeletion_20x_annot.bed) is too small, you cannot see any deletion in the circos images.
Please use file "wen.S1_alldeletion_uniq_annot.bed" to test circosvis program in the example folder.


