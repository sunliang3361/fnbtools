#!/usr/bin/perl 
######## use samtools old version
use warnings; use strict;
use Getopt::Long;
use Time::localtime;
use File::Basename;

#########################  parameters #####################
my $usage = "USAGE:
	$0 
	REQUIRED -c the paired read file 1
	REQUIRED -m the paired read file 2
	

	OPTIONS:
	    -f input annotation file to annotate deletions     
		-h print this help message    


	eg: perl FNBscan.pl -n fnb -c fnb/fnb.mt4_chr1_raw_20x1.bedg -m fnb/fnb.mt4_chr1_mut_20x1.bedg -o fnb/fnb.mt4_chr1_alldeletion_20x.bed -f Mtruncatula_285_Mt4.0v1.gene.gff3
	eg: perl FNBscan.pl -n 20x -c 20x/20x.mt4_sim20x1.bedg -m 20x/20x.mt4_mut300_sim20x1.bedg -o 20x/alldeletion_20x.bed
	eg: perl FNBscan.pl -n 20x_351 -c 20x_351/20x_351.mt4_sim20x1.bedg -m 20x_351/20x_351.mt4_mut351_sim20x1.bedg -o 20x_351/alldeletion_20x_351.bed
	eg: perl FNBscan.pl -n 10x_351 -c 10x_351/10x_351.mt4_sim10x1.bedg -m 10x_351/10x_351.mt4_mut351_sim10x1.bedg -o 10x_351/alldeletion_10x_351.bed
	";

die "$usage\n" if (@ARGV == 0);

my $dir = "/usr/local/fnbtools/";
my $proj   = "fnb"; 
my @cfile;
my @mfile;
my $gff = "";
my $outfile;
my $help = '';
my $cmd;


GetOptions(
    '-c=s@{1,}' => \@cfile,
	'-m=s@{1,}'=> \@mfile,
	'-n=s' =>\$proj,
	'-f=s' =>\$gff,
	'-o=s' =>\$outfile,
	'-h' => \$help
);
# "h"
die "$usage\n" if ($help);

my $result_dirc = $proj; 
my $outfilebase = basename $outfile;
(my $filebase = $outfilebase) =~ s/\.[^.]+$//;


######################### Step1: call VarDiff to identify SNP and small indels##############
my $cfile_vd = process_input(\@cfile,"vcf");
my $mfile_vd = process_input(\@mfile,"vcf");
$cmd = "python ${dir}VarDiff.py -c $cfile_vd -m $mfile_vd -o $result_dirc/$proj.$filebase.unique.vcf";
process_cmd($cmd);

##########################Step2: call DelDiff to identify deletions and merge the deletions identified from clipped and cross reads###############
my $cfiles = join(" ",@cfile);
my $mfiles = join(" ",@mfile);
$cmd = "python ${dir}DelDiff.py -c $cfiles -m $mfiles -f $result_dirc/$proj.$filebase.unique.vcf.del -o $outfile ";
process_cmd($cmd);


##########################Step3: call VarAnnot to annotate deletions (optional)###############
if ($gff){  #we will annotate deletion files
	#annotate all deletions
	(my $ofilebase = $outfile) =~s/\.[^.]+$//;
	my $ofile_del = $ofilebase."_annot.bed";
	#annotate all SNPs
	#my $ofile_snp =
	$cmd = "python VarAnnot.py -i $outfile -t b -f $gff -o $ofile_del";
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

