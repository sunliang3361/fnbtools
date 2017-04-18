#!/usr/bin/perl 
######## use samtools old version
use warnings; use strict;
use Getopt::Long;
#use Statistics::Basic qw(:all);
use List::Util qw[min max];

use Time::localtime;
#use File::Which;

#########################  parameters #####################
my $usage = "USAGE:
	$0
	REQUIRED -s the sam file
	REQUIRED -n the project name
	
		-l library insertion size
		-c <Num,Num,Num> cpu number for 'BWA mem', 'samtools view'  and 'samtools sort', [defualt 8,2,2]
	         
		-h print this help message    

	
		eg: perl  $0 -g genome.fa   -1 reads.fq1 -2 reads.fq2 
		eg: perl FNBscan5 -g mt4.fa -1 ngs-0r7txc096a_S1_L001_R1_001.fastq, ngs-22r0dky60y_S2_L001_R1_001.fastq -2 ngs-0r7txc096a_S1_L001_R2_001.fastq ngs-22r0dky60y_S2_L001_R2_001.fastq
	
	";

die "$usage\n" if (@ARGV == 0);


my $proj   = "fnb";
my $sam;
my $lib_size;
#my $tmp_dir = "$proj.".time();
my $outfile = $proj."infoReads.loc.lst"; #">infor_reads.sam";
my $tmp_dir = "result";
my $help = '';

GetOptions(
	'-l=i' => \$lib_size,
    '-s=s' => \$sam,
	'-n=s' => \$proj,
	'-h' => \$help
);
# "h"
die "$usage\n" if ($help);




########################extract informative reads from sam file##################################
my %reads;
my %dist; #the distance between two paired end reads
my $last;
my $clip_flag = 0;
my $rpos_flag;
my @tlen_numb;
my $tlen_mean;
my $tlen_std;

# #get the mean, standard dev
# open SAM, $sam or die $!;

# while(<SAM>){
	# chomp;
	# next if (/^@/);
	# my ($id, $flag, $chr,$pos,$mapq, $cigar,$rnext,$rpos, $tlen,$seq, $q) = (split /\t/)[0..10];
	# push @tlen_numb, $tlen;
# }

# $tlen_mean = mean(@tlen_numb);
# $tlen_std = stddev(@tlen_numb);


# close(SAM);
open SAM, $sam or die $!;

open OUT, ">infoReads.sam" or die $!;
open LOC, ">$outfile" or die $!;

while(<SAM>){
	chomp;
	next if (/^@/);
	my ($id, $flag, $chr,$pos,$mapq, $cigar,$rnext,$rpos, $tlen,$seq, $q) = (split /\t/)[0..10];
	
	#get all clipped reads and cross reads from sam file
	#output all informative reads to file (with reads ID, chromosome, direction, breakpoint, CRR(cross reads)/CLR(ClippedReads ), Mapping quality)
	#deletion dash line edit from original samtools???????????????????????????????????
	#cluster all informative reads to a deletion files (chromosome, breakpoint,endpoint, deletionlength, SR=crossReads, clipped reads) for each sample
	#confirm deletion by FNBtool and fitler out deletons from control sample.
	if ($last and $id ne $last and %reads){

		#print soft clipped reads 2, will ignore the distance abnormality.
		if ($clip_flag){
			print OUT $reads{$clip_flag}."\n";
		
		}
		#print cross reads
		elsif(defined $dist{1} and defined $dist{2} ){
			if ($dist{1}>= (2*$lib_size) and $dist{2}>= (2*$lib_size)){    #cross read distance larger than 2 library size
				print OUT $reads{1}."\n";
				print OUT $reads{2}."\n";
				my ($crr_id1,$crr_chr1,$crr_pos1) = (split /\t/, $reads{1})[0,2,3];
				my ($crr_id2,$crr_chr2,$crr_pos2) = (split /\t/, $reads{2})[0,2,3];
				my $crr_pos = min($crr_pos1,$crr_pos2);
				my $crr_loc = $crr_pos + ($lib_size/2);
				print LOC "$crr_id1\t$crr_chr1\t$crr_loc\tCRR\n";
			}
		}
		
		undef %reads;
		undef %dist;
		$clip_flag = 0;
		$rpos_flag = 0;
	}
	
	$last = $id;
	$flag = unpack("B32", pack("N",$flag));
	my @bw = (split //, reverse($flag))[0..11]; #bitwise data for flag 
	next if ($bw[2] and $bw[3]);
	my $dir = $bw[6]?1:2;
	#print clipped reads 1 read 1 if read 1,1,2 pattern.  print clipped read2, read2 if read 1, 2,2 pattern.
	if (exists $reads{$dir} and $rpos eq $rpos_flag and $bw[11]){
		print OUT $reads{$dir}."\n";
		print OUT $_."\n";
		$clip_flag = $dir==1?2:1; #if reads1,read1,read2 is the order in sam, print reads2 finally. if reads1, read2, read2, print reads1 finally in the the code line95
		#extract breakpoint of soft clipped read1
		my ($clr_id,$clr_chr,$clr_pos,$clr_cigar)= (split /\t/, $reads{$dir})[0,2,3,5];
		if($clr_cigar=~/^(\d+)M\d+[SH]/){
		#################----------sssssss*******>1  2<sssssssssssss ###############star represent soft cliped
			my $clr_loc = $clr_pos + $1;
			print LOC "$clr_id\t$clr_chr\t$clr_loc\tCLR\n";
		}
		else{
		#################sssssss*******>1 ----------  2<sssssssssssss ###############star represent soft cliped
			my ($clr_id,$clr_chr,$clr_pos,$clr_cigar)= (split /\t/, $_)[0,2,3,5];
			if($clr_cigar=~/^(\d+)M\d+[SH]/){
				my $clr_loc = $clr_pos + $1;
				print LOC "$clr_id\t$clr_chr\t$clr_loc\tCLR\n";
			}else{
				next;
				#print "attention:".$_."\n";
			}
		}

		#edit sam file to get the deletion bar???????????????????????????????????
	}
	else{
		$rpos_flag= $rpos;
		$dist{$dir} = abs($tlen);
		$reads{$dir} = $_;
	}
	
	
	

}
close(SAM);
close(LOC);
close(OUT);

####################################clustering the candidate deletion location###################################
##borrow the idea from itis
my %infoRead;
my @clust;
my @clust_lst;
my $window = $lib_size/2; #need to change so that users can change the size by themselves ???????????????????????????????????????
my $clust_flag = 0;
#sort infoReads.loc.lst file according to the breakpoint position
system("sort -k 2,2 -k 3,3n $outfile >$proj.infoReads.loc.sorted.lst");
open DEL, $proj.".infoReads.loc.sorted.lst" or die $!;


while(<DEL>){
	chomp;
	my($id,$chr,$bk_pos,$typ) = split /\t/;  #$bk_pos is the breakpoint position
	my $win; #determin the clusters
	if (!%infoRead){
		readInfoRead($_,$chr,$bk_pos,$typ); #get the clust info
		next;
	}
	#determin the window size for cross reads and clipped reads
	if ($typ=~/CCR/){
		#cross reads
		$win = $window;
	}elsif($typ=~/CLP/){
		#clipped reads
		$win = 50;
	}else{
		$win = $window;
	}
	
	#get a cluster of support reads
	if($chr eq $infoRead{chr} and abs($bk_pos-$infoRead{pos})<=$win){
		readInfoRead($_,$chr,$bk_pos,$typ);
		if(eof(DEL)){
			#collect read info for this cluster
			my $delInfo = analyzeClust(@clust);
			print "Clust:".$clust_flag."\n";
			$clust_flag++;
			push @clust_lst,$delInfo;
			
		}
	
	
	}else{
		my $delInfo = analyzeClust(@clust);
		print "Clust:".$clust_flag."\n";
		$clust_flag++;
		push @clust_lst,$delInfo;
		
		undef(%infoRead);
		undef(@clust);
		
		readInfoRead($_,$chr,$bk_pos,$typ);
		
		if(eof(DEL)){
			#collect read info for this cluster
			my $delInfo = analyzeClust(@clust);
			print "Clust:".$clust_flag."\n";
			$clust_flag++;
			push @clust_lst,$delInfo;
			
		}
	
	}
	
	
	
}

close DEL;

#write clust into bed file
open BED, ">$proj.bed" or die $!;
writeClust(@clust_lst);

close BED;

#subfunctions

sub readInfoRead{
	my ($it,$chr,$bk_pos,$typ) = @_;
	push @clust, $it;
	$infoRead{chr} = $chr;
	$infoRead{pos} = $bk_pos;
	$infoRead{ty} = $typ;

}

sub analyzeClust{
	my @rcd = @_;
	my @clr_site; # soft clipped read site
	my @crr_site; # cross read site which is not the exact site
	my($id,$chr,$bk_pos,$typ);
	
	foreach my $line (@rcd){
		($id,$chr,$bk_pos,$typ)= split /\t/,$line;
		if ($typ=~/CLR/){
			push @clr_site,$bk_pos;
		}elsif($typ =~/CRR/){
			push @crr_site,$bk_pos;
		}
	}
	#identify the exact deletion site
	my $breakpoint;
	if(@clr_site){
		$breakpoint = mode(@clr_site);   #get the mode value for two or less number????????????????????????????????????????????????????????????????
	}elsif(@crr_site){
		$breakpoint = median(@crr_site);
	}else{
		$breakpoint = 0;
	}
	
	#return useful information
	my ($clr_size,$crr_size) = (scalar@clr_site,scalar@crr_site);
	
	return ("$chr\t$breakpoint\tCLR=$clr_size;CCR=$crr_size;");

}

sub writeClust{
	@clust = @_;
	my $n = 1;
	print BED "DEL#\tChr\tBreakpoint_pos\tSuppRead#\n";
	for my $line (@clust){
		print BED $n."\t".$line."\n";
		$n++;
	}
}
##borrow the idea from itis
sub mode{
	my @array = @_;
	my %rank ;
	foreach my $n (@array){
		$rank{$n}++;
	}
	my $mode = (sort{$rank{$a}<=> $rank{$b}} keys %rank)[-1];
	return $mode;
}
##borrow the idea from itis
sub median {
	my @array = @_;
	@array = sort{$a<=>$b}@array;
	my $size = scalar @array;
	my $median = $array[int($size/2)];
	return $median;
}



