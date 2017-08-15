# Installing FNBTools

## Operating System Requirements

FNBTools has been tested on the following Linux distributions:

* Ubuntu 14.04 LTS
* Ubuntu 16.04 LTS
* CentOS 7.3
* Debian 7 "Wheezy"
* Debian 8 "Jessie"

## Install from Source

1. Run: `git clone ---` (or download fnbtools.tar.gz from Releases and run `tar -xzvf fnbtools.tar.gz`)
2. Run: `sudo ./Install.sh`
3. If prompted with "Would you like to configure as much as possible automatically? [yes]", type **yes** and then press Enter. This will automatically configure Perl's CPAN utility so that additional Perl modules can be installed.
4. If prompted with "Would you like me to automatically choose some CPAN mirror sites for you? [yes]", type **yes** and then press Enter.  This will automatically configure Perl's CPAN utility with a mirror site from which it can download Perl modules.

**Optional** - If you would like to enable visualization of the results, continue with steps 5 - 7:

5. Run: `echo 'export PATH=/usr/local/circos/current/bin:$PATH' >> ~/.bashrc`
6. Run: `. ~/.bashrc`
7. Run `circos -modules` to verify that all Perl modules that are required to run Circos for visualizations are installed. If any are reported as missing, you must install them before attempting to visualize FNBTools results.
    * A helper script called `ReinstallCircosPerlModules.sh` has been provided to assist with installing Perl modules that fail to install during the main installation.
    * Run: `sudo ./ReinstallCircosPerlModules.sh` to re-attempt to install missing Perl modules that Circos requires.

## Install in Docker container

1. git clone ---
2. 


# Using FNBTools

## Step 1 - Align all reads to reference genome

### Usage:

`fnbalign -n fnb -g genome.fa -1 control_forward.fq mutant_forward.fq -2 control_reverse.fq mutant_reverse.fq`

### Parameters:

* REQUIRED -g the genome sequence file in fasta format  
* REQUIRED -1 the paired read file 1
* REQUIRED -2 the paired read file 2
* REQUIRED -n the name of your project (output files will be placed in a directory with the name you provide)

### Example:

An example data set is provided with this repository.

Running the following will create a project directory named 'fnb' relative to where you run the command, and will produce example output from Step 1:

`fnbalign  -n fnb -g example/mt4_chr1_2Mb.fa -1 example/mt4_chr1_raw_20x1.fq example/mt4_chr1_mut_20x1.fq -2 example/mt4_chr1_raw_20x2.fq example/mt4_chr1_mut_20x2.fq`


### Output:

#### 1. BAM files

FNBTools produces one BAM file for each sample.

The output is placed in ./fnb/temp.

Note that the ./fnb directory is named after your project, and was specified as the `-n` switch to the `fnbalign` command.

Running Step 1 using the example data set provided would produce the following BAM files:

* ./fnb/temp/**fnb.mt4_chr1_mut_20x1.ref.sort.bam**
* ./fnb/temp/**fnb.mt4_chr1_raw_20x1.ref.sort.bam**

#### 2. BEDG files

FNBTools produces one Gap file for each sample.

The output is placed in ./fnb (i.e. in the directory named after your project)

Running Step 1 using the example data set provided with this repository would produce the following BEDG files:

* ./fnb/**fnb.mt4_chr1_mut_20x1.bedg**
* ./fnb/**fnb.mt4_chr1_raw_20x1.bedg**

#### 3. BED files

FNBTools produces one BED file for the mutant samples, and one BED file for the control samples.

The output is placed in ./fnb (i.e. in the directory named after your project)

The BED files contain all preliminary deletion data.

Running Step 1 using the example data set provided with this repository would produce the following BED files:

* Mutant samples: ./fnb/**fnb.mt4_chr1_mut_20x1.bed**
* Control samples: ./fnb/**fnb.mt4_chr1_raw_20x1.bed**

## Step 2 - Identify deletions in mutant sample(s)

### Usage:

`fnbscan -n fnb -c control.bedg -m mutant.bedg`

### Parameters:

* REQUIRED -n the project name (use the name you chose in Step 1)
* REQUIRED -c the control BED file
* REQUIRED -m the mutant BED file
* REQUIRED -o output file name containing identified deletions

* OPTIONAL:
  * -f input annotation file to annotate deletions     
  * -h print this help message 
	* -c the control bedg file(s),used to filter out homo and heter deletion
	* -x the contorl bedg file(s),used to filter out homo deletions
	* -g input annotation file to annotate deletions     
	* -h print this help message
	* -a print all homo deletions in mutant including the deletions exist in control * sample [0|1, default:0]
	* -r the overlapping rate between gaps and informative deletion at the same * genomic regions [default:0.9]
	* -d the minimal distance between the breakpoint of informative reads and the * start postion of gap [default:20]
	* -b the minimal crossed reads when there is no clipped reads [default:3]
	* -s the minimal small deletion reads [default:2]
	* -f the minimal flanking reads up and downstream of deletions [default:1]

### Example:

1.	No control sample
`fnbscan -n fnb -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -g example/Mtruncatula_285_Mt4.0v1.gene.gff3`

2.	Filter out all homozygous and heterozygous deletions in control samples
`fnbscan -n fnb -c fnb/fnb.mt4_chr1_raw_20x1.bedg -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -g example/Mtruncatula_285_Mt4.0v1.gene.gff3`

3.	Filter out only homozygous deletions in control samples
`fnbscan -n fnb -x fnb/fnb.mt4_chr1_raw_20x1.bedg -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -g example/Mtruncatula_285_Mt4.0v1.gene.gff3`


### Output:

#### 1. BED file

FNBTools produces a single BED file which contains all unique deletions that were identified.

The output is placed in ./fnb (i.e. in the directory named after your project)

Running Step 2 using the example data set provided with this repository, along with the results of running Step 1 would produce the following BED file:

* ./fnb/**fnb.mt4_chr1_alldeletion_20x.bed**

#### BED file structure

* DEL#: deletion number. All deletions will be sorted by the chromosome number and start postions
* Chr: Chromosome number.	
* BreakpointStart: Start position of deletions
* BreakpointEnd: End position of deletions
* DeletionLength: Deletion length
* SuppRead#: CLR represent the clipped reads number; CRR represents discordant reads number; SMD represents small deletion reads number; FLR represents the flanking supporting reads of deletions on the left; FRR represents the flanking supporting read number of deletions on the right 
* GapStarts_position: Start position of gaps	
* GapEnd_position	: End position of gaps
* Del_mutant: 'Yes' represents deletions found in mutant samples	
* Del_control: 'Yes' represents deletions found in control samples; 'No' represents no deletions found in control samples
* Homo_Unique: 'Yes' represents deletions only exist in mutant samples but not in control sample; 'No' represents deletions exist in both mutant and control sample
* Genes (optional): This column will only show genes if deletions cover

#### 2. BED file with annotation (optional)

If you specify an argument to the `-f` switch which provides an annotation file to annotate deletions, a second BED file with the annotations will be created.

The output is placed in ./fnb (i.e. in the directory named after your project)

Running Step 2 using the example data set provided with this repository, along with the results of running Step 1 **with a `-f` switch argument** would produce the following additional BED file:

* ./fnb/fnb.mt4_chr1_alldeletion_20x_annot.bed

## Step3 - Visualize all identified deletions with Circos

### Usage:

`CircosVis [options] -l <bigDeletion.vcf> -o <output file>`

### Parameters:

	REQUIRED -l the large deletion file
	REQUIRED -o the output image file

### Example:

`circosvis  -l example/wen.S1_alldeletion_uniq_annot.bed -o fnb_circos.png`

`circosvis  -l fnb/fnb.mt4_chr1_alldeletion_20x_annot.bed -o fnb_circos.png`

### Output:

The output circos image will be in 'vis' folder with the name parameter 'o' defined

Since our sample deletion file from step 2 
(./fnb/fnb.mt4_chr1_alldeletion_20x_annot.bed) is too small, you cannot see any deletions in the circos image.

Please use file "wen.S1_alldeletion_uniq_annot.bed" in the example directory of this repository to see example output of the circosvis program.