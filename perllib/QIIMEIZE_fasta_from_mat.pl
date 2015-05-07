#! /usr/bin/perl -w

die "Usage: mat fasta > Redicted\n" unless (@ARGV);
($mat, $new) = (@ARGV);
chomp ($new);
chomp ($mat);

$/=">";

open (IN, "<$new" ) or die "Can't open $new\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    ($info, @seqs) = split ("\n", $line);
    $OTU=$info;
    ($sequence)=join ("", @seqs);
    $seqhash{$OTU}=$sequence;
}
close (IN);

$/="\n";
$first=1;
$counthash=();
open (IN, "<$mat" ) or die "Can't open $mat\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    ($OTU, @pieces)=split ("\t", $line);
    if ($first){
	(@headers)=@pieces;
	$first=0;
    } else {
	if ($seqhash{$OTU}){
	    $tot++;
	    $i=0;
	    $j=@pieces;
	    until ($i >=$j){
		if ($pieces[$i]){
		    $k=0;
		    until ($k >=$pieces[$i]){
			#make a new sequence for each count of the mat for this OTU
			$counthash{$headers[$i]}++;
			print ">$headers[$i]_${counthash{$headers[$i]}} $OTU\n$seqhash{$OTU}\n";
			$k++;
		    }
		}
		$i++;
	    }
	} else {
	    die "Missing sequence $OTU $seqhash{$OTU}\n";
	}
    }
    
}
close(IN);

