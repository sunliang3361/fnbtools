#!/usr/bin/perl 

#Copyright (c) 2017 Noble Research Institute, LLC.
#This portion of FNBTools depends on the following open source software:
#	Samtools, Copyright (C) 2008-2014 Genome Research Ltd., The MIT/Expat License
#	Bcftools, Copyright (C) 2008-2014 Genome Research Ltd., The MIT/Expat License
#	Bedtools, Copyright (C) 1989, 1991 Free Software Foundation, Inc., GNU GENERAL PUBLIC LICENSE
#	BWA, Copyright (C) 2007 Free Software Foundation, Inc., GNU GENERAL PUBLIC LICENSE

use warnings; use strict;
use FindBin;
use Getopt::Long;
use File::Basename;
use Time::localtime;
#use File::Which;

#########################  parameters #####################
my $usage = "USAGE:
	$0
	REQUIRED -g the genome sequence file in fasta format  
	REQUIRED -1 the paired read file 1
	REQUIRED -2 the paired read file 2
	REQUIRED -n the the name of you project

		-p <Num Num Num> cpu number for 'BWA mem', 'samtools view'  and 'samtools sort', [defualt 8 2 2]
	    -l <int> the size of library fragment     
		-h print this help message    

	
		eg: perl fnbalign.pl -n fnb -g rice_chr1_200k.fa -1 rice_chr1_20x_mut_1.fq rice_chr1_20x_raw_1.fq -2 rice_chr1_20x_mut_2.fq rice_chr1_20x_raw_2.fq
		eg: perl fnbalign.pl -n fnb -g mt4_chr1_2Mb.fa -1 mt4_chr1_raw_20x1.fq mt4_chr1_mut_20x1.fq -2 mt4_chr1_raw_20x2.fq mt4_chr1_mut_20x2.fq
		eg: perl fnbalign.pl -n 20x -g mt4.fa -1 mt4_mut300_sim20x1.fq mt4_sim20x1.fq -2 mt4_mut300_sim20x2.fq mt4_sim20x2.fq
		eg: perl fnbalign.pl -n 10x_351 -g mt4.fa -1 mt4_mut351_sim10x1.fq mt4_sim10x1.fq -2 mt4_mut351_sim10x2.fq mt4_sim10x2.fq
		eg: perl fnbalign.pl -n 20x_351 -g mt4.fa -1 mt4_mut351_sim20x1.fq mt4_sim20x1.fq -2 mt4_mut351_sim20x2.fq mt4_sim20x2.fq
		eg: perl fnbalign.pl -n 50x_351 -g mt4.fa -1 mt4_mut351_sim50x1.fq mt4_sim50x1.fq -2 mt4_mut351_sim50x2.fq mt4_sim50x2.fq
	";

die "$usage\n" if (@ARGV == 0);

#my $dir = "/usr/local/fnbtools";
my $dir = "$FindBin::Bin";
my $genome;
my $proj   = "fnb";                 
my @rs1_originals;
my @rs2_originals;
my $lib_len = 500;
my @cpu;
#my $cpu = 8;

my $help = '';

GetOptions(
    '-g=s' => \$genome,
    '-1=s@{1,}' => \@rs1_originals,
    '-2=s@{1,}' => \@rs2_originals,
	'-n=s' =>\$proj,
	'-l=i' =>\$lib_len,
    '-p=s@{3}' => \@cpu,
	#'-p=i' =>\$cpu,
	'-h' => \$help
);

die "$usage\n" if ($help);

my($cpu_bwa,$cpu_view,$cpu_sort) = @cpu;
if(!$cpu_bwa || !$cpu_view|| !$cpu_sort){
	$cpu_bwa = 8;
	$cpu_view = 2;
	$cpu_sort = 2;
}
#my $cpu_bwa= my$cpu_view=my$cpu_sort = $cpu;
#print $cpu_bwa, $cpu_view,$cpu_sort."\n";

my $result_dirc = $proj;
my $result_temp = $proj."/temp";             #we are here



if(scalar(@rs1_originals)!= scalar(@rs2_originals)){
	print "Error: read files doesn't match!\n";
	exit 0;
}
if (!@rs1_originals ||!@rs2_originals){
	print "Your input files do not exist!\n";
	exit 0;
}

########################check exist of required tools##################################
#software must be in the path
####check software
my $samtools;
my $bwa;
my $bcftools;
my $bedtools;
my $cmd;

for my $path (split /:/, $ENV{PATH}){
	if( -f "$path/samtools" && -x _){
		$samtools="$path/samtools";		
	}
		if( -f "$path/bwa" && -x _){
		$bwa="$path/bwa";		
	}
		if( -f "$path/bedtools" && -x _){
		$bedtools="$path/bedtools";		
	}
}
die "No Samtools command available\n" unless ( $samtools );
die "No bwa command available\n" unless ( $bwa );
die "No bedtools command available\n" unless ( $bedtools );



##########################################################

if(-e $result_dirc){
	print "using directory: $result_dirc\n";
}else{
	$cmd = "mkdir $result_dirc";
	system($cmd) == 0 or die $!;
}


if(-e $result_temp ){
	print "using temperary directory: $result_temp \n";
	#empty all temporary bam files
	$cmd = "rm $result_temp/*";
	system($cmd) == 0 or die $!;
}else{
	$cmd = "mkdir $result_temp";
	system($cmd) == 0 or die $!;
}

open CMD, ">$result_temp/run.log" or die $!;



###################### Step 1.1 Index genome sequence ########
#bwa index mt4.fa 

$cmd="cp $genome $result_temp/$proj.ref.fa";
process_cmd($cmd);

$cmd = "bwa index $result_temp/$proj.ref.fa";
process_cmd($cmd);			

$cmd = "samtools faidx $result_temp/$proj.ref.fa";
process_cmd($cmd);	




##### Step 1.2 align original reads to reference genome ######
#bwa mem -t 20  mt4.fa ngs-0r7t3d_S1_L001_R1_001.fastq ngs-0r7t3d_S1_L001_R2_001.fastq >aln.sam
#samtools view
# #samtools sort
#samtools index $result_dirc/$proj.all_reads_aln_reference.sort.bam
my $findex=0;
foreach my $rs1 (@rs1_originals){

	my $rs2=$rs2_originals[$findex];
	my $rs2_file = substr((basename $rs2),0,index((basename $rs2),'.'));
	my $rs1_file=substr((basename $rs1),0,index((basename $rs1),'.'));
	print "input1: $rs1, $rs1_file\n";
	print "input2: $rs2, $rs2_file\n";	

	#my $transformtobed_bam ;
	# if($fast =~ /N/i and $bam == 0){
	#$cmd = "bwa mem -T 20 -t $cpu_bwa $result_temp/$proj.ref.fa $rs1 $rs2 |tee $result_temp/$proj.$rs1_file._ref.sam  |samtools view -@ $cpu_view -buS -q 30 - | samtools sort -@ $cpu_sort  - -O bam -o $result_temp/$proj.$rs1_file._ref.sort.bam";
	$cmd = "bwa mem -T 20 -t $cpu_bwa $result_temp/$proj.ref.fa $rs1 $rs2 > $result_temp/$proj.$rs1_file.ref.sam";
	process_cmd($cmd);
	
	#call extractInfoRead here to extract all informative reads (cross reads and clipped reads) and modify the sam file to get the deletion bar for IGV
	$cmd = "perl -I $dir $dir/extractInfoRead.pl -s $result_temp/$proj.$rs1_file.ref.sam -n $result_dirc/$proj.$rs1_file -l $lib_len "; 
	process_cmd($cmd);
	
	# remove the sam file which takes too much space
	$cmd = "rm $result_temp/$proj.$rs1_file.ref.sam";
	#process_cmd($cmd);
	
	#convert sam to bam, sort and index bam file
	$cmd = "samtools view -@ $cpu_view -buS -q 30 $result_dirc/$proj.$rs1_file.fix.sam | samtools sort -@ $cpu_sort  - -O bam -o $result_temp/$proj.$rs1_file.ref.sort.bam";
	process_cmd($cmd);
	
	$cmd = "samtools index $result_temp/$proj.$rs1_file.ref.sort.bam";
	process_cmd($cmd);

	# remove the fixed sam file which takes too much space
	$cmd = "rm $result_dirc/$proj.$rs1_file.fix.sam";
	#process_cmd($cmd);

	##### Step 2 add supporting reads for flanking reads left and right by using the bam file with mapping score > 30 #####
	open(BED,"<$result_dirc/$proj.$rs1_file.bed") or die $!;
	open(BEDfix,">$result_dirc/$proj.$rs1_file.fix.bed") or die $!;
	my $flag = 0;
	while(my $line = <BED>){
		chomp($line);
		if($flag == 0){
			print BEDfix $line."\n";
			$flag ++;
			next;
		}
		my ($del_n, $chr,$breakpoint,$breakpoint_end,$deletion,$suppRead) = split /\t/, $line;
		my $fl = $breakpoint - $lib_len/2; #$lib_len/3; #the start position of left flanking region
		my $fr = $breakpoint_end + $lib_len/2; #$lib_len/3; #the end position of right flanking region
		#???????????????????????????add flanking reads number here??????????????????????????????????????
		open my $fl_depth, "samtools depth -r $chr:$fl-$breakpoint $result_temp/$proj.$rs1_file.ref.sort.bam |" or die $!;
		open my $fr_depth, "samtools depth -r $chr:$breakpoint_end-$fr $result_temp/$proj.$rs1_file.ref.sort.bam |" or die $!;
		my $fl_n = 0;
		my $fr_n = 0;
		
		while(<$fl_depth>){
			chomp;
			my @records = split /\t/,$_;
			$fl_n = $fl_n + $records[2];
		}
		while(<$fr_depth>){
			chomp;
			my @records = split /\t/,$_;
			$fr_n = $fr_n + $records[2];
		}
		my $flr = int($fl_n*2/$lib_len);
		my $frr = int($fr_n*2/$lib_len);
		# my $flr = int($fl_n/20);
		# my $frr = int($fr_n/20);
		print BEDfix $line."FLR=$flr;FRR=$frr\n";
		close $fl_depth;
		close $fr_depth;
	}
	close(BED);
	close(BEDfix);

	# get the gap file
	#$cmd= "genomeCoverageBed -ibam $result_temp/$proj.$rs1_file.ref.sort.bam -bga |";
	$cmd= "bedtools genomecov -ibam $result_temp/$proj.$rs1_file.ref.sort.bam -bga |";
	open (OUT, ">$result_dirc/$proj.$rs1_file.bedg") or die ("cannot open bedg file");
	open IN, $cmd;
	while(my $line =<IN>){
		my @cols = split("\t",$line);
		if ($cols[3] == 0){
			print OUT $line;
		}
	
	}
	close(IN);
	close(OUT);
	

	
	$findex=$findex+1;
}



	
#  subfuctions derived  from trinity package
sub process_cmd {
    my ($cmd) = @_;

	print STDERR &mytime."CMD: $cmd\n";
	print CMD "$cmd\t#".&mytime."\n";
	my $start_time = time();
	my $ret = system($cmd);
	my $end_time = time();

	if ($ret) {
		die "Error, cmd: $cmd died with ret $ret";
	}

	print STDERR "CMD finished (" . ($end_time - $start_time) . " seconds)\n";    

	return;

}


sub mytime() {
  my @mabbr = qw(January February March April May June July August September October November December);
  my @wabbr = qw(Sunday Monday Tuesday Wednesday Thursday Friday Saturday);
  my $sec = localtime->sec() < 10 ? '0' . localtime->sec() : localtime->sec();
  my $min = localtime->min() < 10 ? '0' . localtime->min() : localtime->min();
  my $hour = localtime->hour() < 10 ? '0' . localtime->hour() : localtime->hour();
  my $wday = $wabbr[localtime->wday];
  my $mday = localtime->mday;
  my $mon = $mabbr[localtime->mon];
  my $year = localtime->year() + 1900;
  return "$wday, $mon $mday, $year: $hour:$min:$sec\t";
}
