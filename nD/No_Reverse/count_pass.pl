#! /usr/bin/perl
use strict;

my $sam = $ARGV[0];

$sam =~ /(.+)\.sam/;
my $prefix = $1;
open OUT, ">$prefix.reads.info" or die $!;

my (%MIS_POS, %DEL_POS, %INS_POS, %ERR_POS);

my (%INS,%DEL);

my %REAL_LEN;
print OUT "READ_ID\tMatch\tMismatch\tInsertion\tDeletion\tLen\tStart\n";
open IN, $sam or die $!;
my %hash;
my ($total_match, $total_mismatch, $total_ins, $total_del, $total_len) = (0,0,0,0,0);
while (<IN>){
	next if /^@/;
	chomp;
	my ($qname, $flag, $cigar, $seq) = (split /\s+/, $_)[0,1,5,9];
	#print $qname, "\n";
	my ($match, $mismatch, $ins, $del)=(0,0,0,0);
	next if $cigar eq "*";

#=================
	my ($head_soft, $head_hard, $tail_soft,$tail_hard) = (0,0,0,0);
	if ($cigar=~ /^(\d+)[S|H]/){  #head
		if ($cigar=~ /^(\d+)S/){
			$head_soft = $1;
		}elsif ($cigar=~ /^(\d+)H/){
			$head_hard = $1;
		}
	}
	if ($cigar=~ /(\d+)[S|H]$/){ #tail
		if ($cigar=~ /(\d+)S$/){
		        $tail_soft = $1;
		}elsif ($cigar=~ /(\d+)H$/){
			$tail_hard = $1;
		}
	}
	my $len = (length $seq)-$head_soft -$tail_soft;
	my $real_len = (length $seq) + $head_hard + $tail_hard;
	$REAL_LEN{$qname} = $real_len;
	
	my $beg =0;
	if ($flag =~ /[0|256|2048]/){ # forward mapping
	       	if ($cigar =~ /^(\d+)[S|H]/){
			$beg = $1;	
		}
	}elsif ($flag =~ /[16|272|2064]/){ #reverse mapping
		if ($cigar =~ /(\d+)[S|H]$/){
			$beg = $1;
		}
	}
	
#	my ($head, $tail) = (0,0); 
#	if ($cigar=~ /^(\d+)S/){
 #       	$head = $1;
#	}
#	if ($cigar=~ /(\d+)S$/){
 #       	$tail = $1;
#	}	
#	my $len = (length $seq)-$head-$tail;
#=================
	my $acu_base=0;
	my $win;
	while ($cigar !~ /^$/){
		if ($cigar =~ /^([0-9]+[=XIDSH])/){
			my $cigar_part = $1;
			if ($flag =~ /[0|256|2048]/){ # forward mapping
        			$win = int(($beg + $acu_base)/5)+1;
			}elsif ($flag =~ /[16|272|2064]/){ #reverse mapping
				$acu_base = $len - $acu_base;
				$win = int(($beg + $acu_base)/5)+1;
			}

			if ($cigar_part=~ /^(\d+)=/){
				$match += $1;
				$acu_base += $1;
			}elsif($cigar_part =~ /(\d+)X/){
				$mismatch += $1;
				$acu_base += $1;
				$MIS_POS{$qname}{$win}+=$1;
				$ERR_POS{$qname}{$win}+=$1;
			}elsif ($cigar_part =~ /(\d+)I/){
				$ins +=$1;
				$acu_base += $1;
				$INS{$1}++;
				$INS_POS{$qname}{$win}+=$1;
				$ERR_POS{$qname}{$win}+=$1;
			}elsif ($cigar_part =~ /(\d+)D/){
				$del +=$1;
				$DEL{$1}++;
				$DEL_POS{$qname}{$win}+=$1;
				$ERR_POS{$qname}{$win}+=$1;
			}
			$cigar =~ s/$cigar_part//;
		}else{
			die "Unexpected cigar: $qname\t$cigar\n";
		}
	}
#	print $qname, "\t", $match, "\t",$mismatch, "\t", $ins, "\t", $del, "\n";
#	print $qname, "\t", $match, "\t",$mismatch, "\t", $ins, "\t", $del, "\t", $len, "\n";
	my $error_rate = ($mismatch+$ins+$del)/($mismatch+$ins+$del+$match);
	my $iden_rate = 1-$error_rate;
	if ($iden_rate >=0.999){
		$hash{">=0.999"}++;
	}elsif ($iden_rate >=0.99){
		$hash{"0.99~0.999"}++;
	}elsif ($iden_rate >=0.98){
		$hash{"0.98~0.99"}++;
	}elsif ($iden_rate >=0.97){
		$hash{"0.97~0.98"}++;
	}elsif ($iden_rate >=0.96){
		$hash{"0.96~0.97"}++;
	}elsif ($iden_rate >=0.95){
		$hash{"0.95~0.96"}++;
	}elsif ($iden_rate >=0.94){
		$hash{"0.94~0.95"}++;
	}elsif ($iden_rate >=0.93){
		$hash{"0.93~0.94"}++;
	}elsif ($iden_rate >=0.92){
		$hash{"0.92~0.93"}++;
	}elsif ($iden_rate >=0.91){
		$hash{"0.91~0.92"}++;
	}elsif ($iden_rate >=0.90){
		$hash{"0.90~0.91"}++;
	}else{
		$hash{"<0.90"}++;	
	}
	
	$total_match += $match;
	$total_mismatch += $mismatch;
	$total_ins += $ins;
	$total_del += $del;
	$total_len += $match+$mismatch+$ins+$del;
	print OUT $qname, "\t", $match, "\t",$mismatch, "\t", $ins, "\t", $del, "\t", $len, "\t", $beg, "\t", $iden_rate,"\n";
}
close OUT;
close IN;

open STAT, ">$prefix.total.stat" or die $!;
my $total_mismatch_rate = $total_mismatch/$total_len;
my $total_ins_rate = $total_ins/$total_len;
my $total_del_rate = $total_del/$total_len;
my $total_iden_rate = 1-$total_mismatch_rate-$total_ins_rate-$total_del_rate;
print STAT "total_mismatch_rate: $total_mismatch_rate\n";
print STAT "total_ins_rate: $total_ins_rate\n";
print STAT "total_del_rate: $total_del_rate\n";
print STAT "total_iden_rate: $total_iden_rate\n";

foreach my $range (sort keys %hash){
	print STAT $range, "\t", $hash{$range},"\n";
}

close STAT;

open INDEL, ">$prefix.indel.stat" or die $!;
print INDEL "Length\tInsertion\tDeletion\n";
foreach my $len (sort {$a<=>$b} keys %INS){
	print INDEL $len, "\t", $INS{$len}, "\t", $DEL{$len}, "\n";
}
close INDEL;


open MISPOS, ">$prefix.reads.mismatch.detail" or die $!;
open INSPOS, ">$prefix.reads.insert.detail" or die $!;
open DELPOS, ">$prefix.reads.delete.detail" or die $!;
open ERRPOS, ">$prefix.reads.error.detail" or die $!;

print ERRPOS "READ_ID\t";
for my $win (1..999){
	print ERRPOS $win*5, "\t";
}
print ERRPOS "5000\n";

foreach my $qname (keys %REAL_LEN){
	print ERRPOS $qname,"\t";
	for my $win (1..999){
		my $error_num;
		if ($ERR_POS{$qname}{$win}){
			$error_num = $ERR_POS{$qname}{$win};
		}else{
			if ($win*5 < $REAL_LEN{$qname}){
				$error_num =0;
			}else{
				$error_num = "NULL";
			}
		}
		print ERRPOS $error_num, "\t";
	}
	print ERRPOS $ERR_POS{$qname}{1000},"\n";
}
close MISPOS;
close INSPOS;
close DELPOS;
close ERRPOS;
