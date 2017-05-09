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

	
		eg: perl  $0 extractInfoRead.pl -s infoReads.sam -n loc.lst.test -l 500
		
	";

die "$usage\n" if (@ARGV == 0);


my $proj   = "fnb";
my $sam;
my $lib_size;
#my $tmp_dir = "$proj.".time();

#my $tmp_dir = "result";
my $help = '';

GetOptions(
	'-l=i' => \$lib_size,
    '-s=s' => \$sam,
	'-n=s' => \$proj,
	'-h' => \$help
);
# "h"
die "$usage\n" if ($help);
my $outfile = $proj.".infoReads.loc.lst"; #
my $outfile_smDel = $proj.".smDel.loc.lst"; #
my $sam_m = $proj.".fix.sam";


########################extract informative reads from sam file##################################
my %reads;
my %dist; #the distance between two paired end reads
my $last;
my $clip_flag = 0;
my %rpos_flag;
my %rnext_flag;
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
open SAMm, ">$sam_m" or die $!;

 
#open OUT, ">infoReads.sam" or die $!;
open LOC, ">$outfile" or die $!;
#open LOC_smDel, ">$outfile_smDel" or die $1;

while(<SAM>){
	chomp;
	if (/^@/){
		print SAMm $_."\n";
		next;
	}
	my ($id, $flag, $chr,$pos,$mapq, $cigar,$rnext,$rpos, $tlen,$seq, $q) = (split /\t/)[0..10];
	#get small deletion from cigar
	my %sm_del = cigar2del($cigar);
	if(%sm_del and $mapq >= 30){
		foreach my $key (keys %sm_del){
			my $sm_pos_st = ($pos + $key - 1);
			my $sm_pos_ed = ($pos + $key -1 + $sm_del{$key});
			print LOC "$id\t$chr\t$sm_pos_st\t$sm_pos_ed\t$sm_del{$key}\tSMD\n";
			#print LOC_smDel "$id\t$chr\t$sm_pos_st\t$sm_pos_ed\t$sm_del{$key}\tSMD\n";
			#print "$id\t$chr\t$sm_pos_st\t$sm_pos_ed\t$sm_del{$key}\tSMD\n";
		}

	}
	
	#get all clipped reads and cross reads from sam file
	#output all informative reads to file (with reads ID, chromosome, direction, breakpoint, CRR(cross reads)/CLR(ClippedReads ), Mapping quality)
	#deletion dash line edit from original samtools???????????????????????????????????
	#cluster all informative reads to a deletion files (chromosome, breakpoint,endpoint, deletionlength, SR=crossReads, clipped reads) for each sample
	#confirm deletion by FNBtool and fitler out deletons from control sample.
	if ($last and $id ne $last and %reads){

		#print soft clipped read 2 if read1 read1 read2. elif read1 read2 read2 print read 1
		if ($clip_flag && exists $reads{$clip_flag}){
			print SAMm $reads{$clip_flag}."\n";
		}
		#print cross reads
		elsif(defined $dist{1} and defined $dist{2} ){
			if ($dist{1}>= (2*$lib_size) and $dist{2}>= (2*$lib_size)){    #cross read distance larger than 2 library size
				#print OUT $reads{1}."\n";
				#print OUT $reads{2}."\n";
				print SAMm $reads{1}."\n";
				print SAMm $reads{2}."\n";
				my ($crr_id1,$crr_chr1,$crr_pos1) = (split /\t/, $reads{1})[0,2,3];
				my ($crr_id2,$crr_chr2,$crr_pos2) = (split /\t/, $reads{2})[0,2,3];
				my $crr_pos = min($crr_pos1,$crr_pos2);
				my $crr_loc = $crr_pos + ($lib_size/2);
				my $crr_loc_end = max($crr_pos1,$crr_pos2);
				my $crr_d = abs($crr_pos1-$crr_pos2) -350;
				
				if($crr_d <= 1500000){
					print LOC "$crr_id1\t$crr_chr1\t$crr_loc\t$crr_loc_end\t$crr_d\tCRR\n";
				}
				
				
			}else{ #not crossed paired end reads
				print SAMm $reads{1}."\n";
				print SAMm $reads{2}."\n";
			}
		}
		# there is no informative reads?????? Just print it out to SAMm
		else{ #single reads here
			if (defined $dist{1}){
				#print "checking this reads:".$reads{1}."\n";
				print SAMm $reads{1}."\n";
			}else{
				print SAMm $reads{2}."\n";
			}

		}
		
		undef %reads;
		undef %dist;
		$clip_flag = 0;
		undef %rpos_flag;
		undef %rnext_flag;
	}
	
	$last = $id;
	$flag = unpack("B32", pack("N",$flag));
	my @bw = (split //, reverse($flag))[0..11]; #bitwise data for flag 
	if ($bw[2] and $bw[3]){
		print SAMm $_."\n";
		next;
	}
	
	# if ($mapq == 0){
		# print SAMm $_."\n";
		# next;
	# }
	
	my $dir = $bw[6]?1:2;
	#print clipped reads 1 read 1 if read 1,1,2 pattern.  print clipped read2, read2 if read 1, 2,2 pattern.
	if (exists $reads{$dir}){
		#print OUT $reads{$dir}."\n";
		#print OUT $_."\n";
		$clip_flag = $dir==1?2:1; #if reads1,read1,read2 is the order in sam, print reads2 finally. if reads1, read2, read2, print reads1 finally in the the code line95
		if ($rnext_flag{$dir} ne '=' || $rnext ne '=' || $rpos != $rpos_flag{$dir} || !$bw[11]) {
			print SAMm $reads{$dir}."\n";
			print SAMm $_."\n";
			if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
			undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
			undef %dist;
			undef %rpos_flag;
			undef %rnext_flag;
			next;
		}
		
		my $clr_cigar_new;
		my $qual_new;
		my $seq_new;
		my $d;
		#extract breakpoint of soft clipped read1
		my ($clr_id,$clr_flag,$clr_chr,$clr_pos,$clr_cigar,$clr_rpos,$seq, $qual)= (split /\t/, $reads{$dir})[0,1,2,3,5,7,9,10];
		my ($clr_id2,$clr_flag2,$clr_chr2,$clr_pos2,$clr_cigar2,$clr_rpos2,$seq2,$qual2)= (split /\t/, $_)[0,1,2,3,5,7,9,10];
		
		if($clr_cigar=~/^(\d+)M\d+[SH]$/){
		################# Read1: 126M24S  READ1:125H25M READ2:150M *OR* READ1:150M Read2: 126M24S  READ2:125H25M 
			my $m1 = $1;
			my $s1 = substr $seq, 0, $m1;   # match sequences
			my $clr_loc = $clr_pos + $m1 -1;
			my $q1 = substr $qual, 0, $m1;

			#change the cigar of clipped reads
			if($clr_cigar2=~/^(\d+)([SH])(\d+)M$/g){
				my $ss_n = $1;    #soft seqs number
				if ($2 eq 'S'){
					my $s2 = substr $seq2, $ss_n,$3;
					my $q2 = substr $qual2,$ss_n,$3;
					$seq_new = $s1.$s2;
					$qual_new = $q1.$q2;
				}else{
					my $s2 = substr $seq2, 0, $3;
					my $q2 = substr $qual2, 0, $3; 
					$seq_new = $s1.$s2;
					$qual_new = $q1.$q2;
				}
				my $m2 = $3;
				#$d = abs($clr_pos2-$clr_pos) - $m1;
				if(abs($ss_n - $m1)>10){ #e.g. 50M100S  40S101M, 50-40<10, good
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
				}
				
				if ($clr_pos2-$clr_pos<=0){ #insertion
					#print "Read_MS\n";
					#print $clr_id2."\n";
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
					#die();
				}
				elsif($clr_pos2-$clr_pos - $m1 < 0){
					#print "Read_MS_m\n";
					#print $clr_id2."\n";
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";	
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
				}
				elsif($clr_pos2-$clr_pos - $m1 > 1500000){
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";	
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
				}
				else{
					$d = $clr_pos2-$clr_pos- $m1;
				}
				$clr_cigar_new = $m1."M".$d."D".$m2."M";
				print SAMm "$clr_id\t$clr_flag\t$clr_chr\t$clr_pos\t60\t$clr_cigar_new\t=\t$clr_rpos\t500\t$seq_new\t$qual_new\n";
				#print informative reads
				#print LOC "$clr_id\t$clr_chr\t$clr_loc\tCLR\n";
				my $clr_loc_end = $clr_loc + $d;
				print LOC "$clr_id\t$clr_chr\t$clr_loc\t$clr_loc_end\t$d\tCLR\n"; ######?????????????????????????????????????????check the same chromose for clipped reads
			}else{
				print SAMm $reads{$dir}."\n";
				print SAMm $_."\n";
				if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
				undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
				undef %dist;
				undef %rpos_flag;
				undef %rnext_flag;
				next;
			}
			
			

			
		}
		elsif($clr_cigar2=~/^(\d+)M\d+[SH]$/){
		#################READ1:37S113M  READ1:37M113H READ2:150M *OR* READ1:150M READ2:37S113M  READ2:37M113H 
			my $m1 = $1;
			my $clr_loc2= $clr_pos2 + $m1 - 1;
			#print LOC "$clr_id2\t$clr_chr2\t$clr_loc2\tCLR\n";
			#Change the cigar of clipped reads
			my $s1 = substr $seq2, 0, $m1;
			my $q1 = substr $qual2, 0, $m1;
		
			if ($clr_cigar=~/^(\d+)([SH])(\d+)M$/g){
				my $ss_n = $1;    #soft seqs number
				if ($2 eq 'S'){
					my $s2 = substr $seq,$ss_n,$3;
					my $q2 = substr $qual,$ss_n,$3;
					$seq_new = $s1.$s2;
					$qual_new = $q1.$q2;
				}else{
					my $s2 = substr $seq, 0, $3;
					my $q2 = substr $seq, 0, $3;
					$seq_new = $s1.$s2;
					$qual_new = $q1.$q2;
				}
				my $m2 = $3;
				
				if(abs($ss_n - $m1)>10){ #e.g. 50M100S  40S101M, 50-40<10, good
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
				}
				#$d = abs($clr_pos - $clr_pos2) - $m1;
				if ($clr_pos-$clr_pos2<=0){  #insertion
					#print "Read_SM\n";
					#print $clr_id2."\n";
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
					#die();
				}
				elsif($clr_pos - $clr_pos2 - $m1 < 0){
					#print "Read_SMm\n";
					#print $clr_id2."\n";
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
				}elsif($clr_pos - $clr_pos2 - $m1>1500000){
					print SAMm $reads{$dir}."\n";
					print SAMm $_."\n";
					if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
					undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
					undef %dist;
					undef %rpos_flag;
					undef %rnext_flag;
					next;
				}
				else{
					$d = $clr_pos - $clr_pos2 - $m1;
				}
				$clr_cigar_new = $m1."M".$d."D".$m2."M";
				print SAMm "$clr_id2\t$clr_flag2\t$clr_chr2\t$clr_pos2\t60\t$clr_cigar_new\t=\t$clr_rpos2\t500\t$seq_new\t$qual_new\n";
				#print LOC "$clr_id\t$clr_chr\t$clr_loc\tCLR\n";
				
				my $clr_loc_end2 = $clr_loc2 + $d;
				#print LOC "$clr_id\t$clr_chr\t$clr_loc\t$clr_loc_end\t$d\tCLR\n";
				print LOC "$clr_id2\t$clr_chr2\t$clr_loc2\t$clr_loc_end2\t$d\tCLR\n";
			}else{
				print SAMm $reads{$dir}."\n";
				print SAMm $_."\n";
				if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
				undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
				undef %dist;
				undef %rpos_flag;
				undef %rnext_flag;
				next;
			}
		}
		else{
			print SAMm $reads{$dir}."\n";
			print SAMm $_."\n";
			if(exists $reads{$clip_flag}){print SAMm $reads{$clip_flag}."\n"; }
			undef %reads; #in case of situtation in sam file, read1 read1 read1 read1 read 2 read2
			undef %dist;
			undef %rpos_flag;
			undef %rnext_flag;
			next;
		
		}

		#edit sam file to get the deletion bar???????????????????????????????????
	}
	else{
		$rpos_flag{$dir} = $rpos;
		$rnext_flag{$dir} = $rnext;
		$dist{$dir} = abs($tlen); 
		$reads{$dir} = $_;
		
	}
	
	
	

}
close(SAM);
close(LOC);
#close(LOC_smDel);
close(SAMm);
#close(OUT);

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
	my($id,$chr,$bk_pos,$bk_end,$del,$typ) = split /\t/;  #$bk_pos is the breakpoint position
	my $win; #determin the clusters
	if (!%infoRead){
		readInfoRead($_,$chr,$bk_pos,$bk_end,$del,$typ); #get the clust info
		next;
	}
	#determin the window size for cross reads and clipped reads
	if ($typ=~/CRR/){
		#cross reads
		$win = $window;
	}elsif($typ=~/CLR/){
		#clipped reads
		$win = 50;
	}elsif($typ=~/SMD/){
		$win = 5;
	}
	else{
		$win = $window;
	}
	#if last reads is small deletion, then the next reads should use a small window even if the next read is clipped reads.
	if($infoRead{ty}=~/SMD/){
		$win = 5;
	}
	
	#get a cluster of support reads
	if($chr eq $infoRead{chr} and abs($bk_pos-$infoRead{pos})<=$win){
		readInfoRead($_,$chr,$bk_pos,$bk_end,$del,$typ);
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
		
		readInfoRead($_,$chr,$bk_pos,$bk_end,$del,$typ);
		
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
sub cigar2del{
	my $cigar = $_[0];
	my %sm_del;
	my $p_len = 0;
	foreach my $p ($cigar=~ /\d+[MIDNSHP=X]{1}/g){
		my $char = chop $p;
		my $num = $p;
		if ($char eq "D"){
			$sm_del{$p_len} = $num;
		}
		if ($char eq "D" or $char eq "M" ){
			$p_len = $p_len + $num;
		}
	}
	return %sm_del;

}

sub readInfoRead{
	my ($it,$chr,$bk_pos,$bk_end,$del,$typ) = @_;
	push @clust, $it;
	$infoRead{chr} = $chr;
	$infoRead{pos} = $bk_pos;
	$infoRead{ty} = $typ;

}

sub analyzeClust{
	my @rcd = @_;
	my @clr_site; # soft clipped read site
	my @crr_site; # cross read site which is not the exact site
	my @smd_site; # small deletion site
	my @clr_site_end;
	my @crr_site_end;
	my @smd_site_end;
	my @clr_site_del;
	my @crr_site_del;
	my @smd_site_del; #deletion length
	my($id,$chr,$bk_pos,$bk_end,$del,$typ);

	
	foreach my $line (@rcd){
		($id,$chr,$bk_pos,$bk_end,$del,$typ)= split /\t/,$line;
		if ($typ=~/CLR/){
			push @clr_site,$bk_pos;
			push @clr_site_end,$bk_end;
			push @clr_site_del,$del;
		}elsif($typ =~/CRR/){
			push @crr_site,$bk_pos;
			push @crr_site_end,$bk_end;
			push @crr_site_del,$del;
		}elsif($typ =~/SMD/){
			push @smd_site,$bk_pos;
			push @smd_site_end,$bk_end;
			push @smd_site_del,$del;
		}
	}

	#return useful information
	my ($clr_size,$crr_size,$smd_size) = (scalar@clr_site,scalar@crr_site,scalar@smd_site);
	
	#identify the exact deletion site
	my $breakpoint;
	my $breakpoint_end;
	my $deletion;
	
	if(@smd_site and @clr_site){
	  if ($smd_size >= $clr_size){
		$breakpoint = mode(@smd_site);
		$breakpoint_end = mode(@smd_site_end);
		$deletion = mode(@smd_site_del);
	  }else{
		$breakpoint = mode(@clr_site);
		$breakpoint_end = mode(@clr_site_end);
		$deletion = mode(@clr_site_del);
	  }
	}
	elsif(@smd_site){
		$breakpoint = mode(@smd_site);
		$breakpoint_end = mode(@smd_site_end);
		$deletion = mode(@smd_site_del);
	}
	elsif(@clr_site){
		$breakpoint = mode(@clr_site);   #get the mode value for two or less number????????????????????????????????????????????????????????????????
		$breakpoint_end = mode(@clr_site_end);
		$deletion = mode(@clr_site_del);
	}
	elsif(@crr_site){
		$breakpoint = median(@crr_site);
		$breakpoint_end = median(@crr_site_end);
		$deletion = median(@crr_site_del);
	}else{
		$breakpoint = 0;
		$breakpoint_end = 0;
		$deletion = 0;
	}
	
	
	return ("$chr\t$breakpoint\t$breakpoint_end\t$deletion\tCLR=$clr_size;CRR=$crr_size;SMD=$smd_size;");

}

sub writeClust{
	@clust = @_;
	my $n = 1;
	print BED "DEL#\tChr\tBreakpointStart\tBreakpointEnd\tDeletionLength\tSuppRead#\n";
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



