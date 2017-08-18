#!/usr/bin/perl 

# Noble Research Institute, LLC License
 
# Copyright (c) 2017 Noble Research Institute, LLC
 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE. UNLESS OTHERWISE PROHIBITED BY LAW, RECIPIENT AGREES TO INDEMNIFY, 
# DEFEND WITH COUNSEL ACCEPTABLE TO PROVIDER, AND HOLD HARMLESS PROVIDER AND 
# PROVIDER’S TRUSTEES, OFFICERS, AGENTS, AND EMPLOYEES FROM ANY AND ALL CLAIMS.

######## use samtools old version
use warnings; use strict;
use Getopt::Long;
use FindBin;
use Time::localtime;
use File::Basename;

#########################  parameters #####################
my $usage = "USAGE:
	$0
	REQUIRED -n the project name
	REQUIRED -m the mutant bedg file(only one pooled sample allowed)
	REQUIRED -o output file name containing identified deletions

	OPTIONS:
	-c the control bedg file(s),used to filter out homo and heter deletion
	-x the contorl bedg file(s),used to filter out homo deletions
	-g input annotation file to annotate deletions     
	-h print this help message
	-a print all homo deletions in mutant including the deletions exist in control sample [0|1, default:0]
	-r the overlapping rate between gaps and informative deletion at the same genomic regions [default:0.9]
	-d the minimal distance between the breakpoint of informative reads and the start postion of gap [default:20]
	-b the minimal crossed reads when there is no clipped reads [default:3]
	-s the minimal small deletion reads [default:2]
	-f the minimal flanking reads up and downstream of deletions [default:1]

	eg: perl fnbscan.pl -n fnb -c fnb/fnb.mt4_chr1_raw_20x1.bedg -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -g Mtruncatula_285_Mt4.0v1.gene.gff3
	eg: perl fnbscan.pl -n 20x -c 20x/20x.mt4_sim20x1.bedg -m 20x/20x.mt4_mut300_sim20x1.bedg -o 20x/alldeletion_20x.bed
	eg: perl fnbscan.pl -n 20x_351 -c 20x_351/20x_351.mt4_sim20x1.bedg -m 20x_351/20x_351.mt4_mut351_sim20x1.bedg -o 20x_351/alldeletion_20x_351.bed
	eg: perl fnbscan.pl -n 10x_351 -c 10x_351/10x_351.mt4_sim10x1.bedg -m 10x_351/10x_351.mt4_mut351_sim10x1.bedg -o 10x_351/alldeletion_10x_351.bed
	";

die "$usage\n" if (@ARGV == 0);

#my $dir = "/usr/local/fnbtools";
my $dir = "$FindBin::Bin";
my $proj   = "fnb"; 
my @cfile;
my @xfile;
my @mfile;
my $gff = "";
my $outfile;
my $help = '';
my $orate = 0.9;
my $minDiff = 20;
my $minCRR = 3;
my $minSMD = 2;
my $minFR = 1;
my $allHomo = 0;
my $cmd;


GetOptions(
    '-c:s@{,}' => \@cfile,
	'-x:s@{,}' => \@xfile,
	'-m=s@{1,}'=> \@mfile,
	'-n=s' =>\$proj,
	'-g=s' =>\$gff,
	'-r=f' =>\$orate,
	'-d=i' =>\$minDiff,
	'-b=i' =>\$minCRR,
	'-s=i' =>\$minSMD,
	'-f=i' =>\$minFR,
	'-a=i' =>\$allHomo,
	'-o=s' =>\$outfile,
	'-h' => \$help
);
# "h"
die "$usage\n" if ($help);

my $result_dirc = $proj; 
my $outfilebase = basename $outfile;
(my $filebase = $outfilebase) =~ s/\.[^.]+$//;


######################### Step1: call VarDiff to identify SNP and small indels##############
# my $cfile_vd = process_input(\@cfile,"vcf");
# my $mfile_vd = process_input(\@mfile,"vcf");
# $cmd = "python VarDiff.py -c $cfile_vd -m $mfile_vd -o $result_dirc/$proj.$filebase";
# process_cmd($cmd);

##########################Step2: call DelDiff to identify deletions and merge the deletions identified from clipped and cross reads###############
if (!@mfile){
	print "Your your mutant files do not exist!\n";
	exit 0;
}


if(@cfile){
	my $mfiles = join(" ",@mfile);
	my $cfiles = join(" ",@cfile);
	$cmd = "python $dir/DelDiff.py -c $cfiles -m $mfiles -f $minFR -r $orate -d $minDiff -b $minCRR -s $minSMD -a $allHomo -o $outfile ";
	process_cmd($cmd);
}
elsif(@xfile){
	my $mfiles = join(" ",@mfile);
	my $xfiles = join(" ",@xfile);
	$cmd = "python $dir/DelDiff.py -x $xfiles -m $mfiles -f $minFR -r $orate -d $minDiff -b $minCRR -s $minSMD -a $allHomo -o $outfile ";
	process_cmd($cmd);
}
else{
	my $mfiles = join(" ",@mfile);
	$cmd = "python $dir/DelDiff.py -m $mfiles -f $minFR -r $orate -d $minDiff -b $minCRR -s $minSMD -a $allHomo -o $outfile ";
	process_cmd($cmd);
}




##########################Step3: call VarAnnot to annotate deletions (optional)###############
if ($gff){  #we will annotate deletion files
	#annotate all deletions
	(my $ofilebase = $outfile) =~s/\.[^.]+$//;
	my $ofile_del = $ofilebase."_annot.bed";
	#annotate all SNPs
	#my $ofile_snp =
	$cmd = "python $dir/VarAnnot.py -i $outfile -t b -f $gff -o $ofile_del";
	process_cmd($cmd);
}


#	subfunctions
sub process_input{
	my($file,$suffix) = @_;
	my $newfile="";
	foreach my $f (@$file){
		#print $f."\n";
		if ($f=~/(.*)\.bedg$/){
			$newfile = $newfile." ".$1.".".$suffix;
			#print $newfile."\n";
		}else{
			die("Error: your input file type is not correct!")
		}
	}
	return $newfile;
}
	
#	subfuctions derived  from trinity package
sub process_cmd {
    my ($cmd) = @_;

	print STDERR &mytime."CMD: $cmd\n";
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

