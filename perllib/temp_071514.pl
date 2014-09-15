#! /usr/bin/perl -w

die "Usage: q.derep.fst q.index output\n" unless (@ARGV);
($derep, $index, $output) = (@ARGV);
chomp ($derep, $index, $output);

$/=">";
open (OUT, ">${output}.fa") or die "Can't open ${output}.fa\n";
open (IN, "<$derep" ) or die "Can't open $derep\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split ("\n", $line);
    ($head)=shift(@pieces);
    ($seqno, $counts)=split (";", $head);
    ($sequence)=join ("", @pieces);
    if ($hash{$seqno}){
	die "Not unique $sequence $hash{$seqno}\n";
    } else {
	$hash{$seqno}=$sequence;
	$countshash{$seqno}=$counts;
	print OUT ">$seqno\n$sequence\n";
    }
}

close (IN);
close (OUT);


$/="\n";
open (IN, "<$index" ) or die "Can't open $index\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    ($lib, $seq, $counts)=split ("\t", $line);
    if ($hash{$seq}){
	$indexhash{$seq}{$lib}=$counts;
	$libhash{$lib}++;
    } else {
	die "Missing value in derep $seq\n";
    }
}
close (IN);

open (MAT, ">${output}.f0.mat") or die "Can't open ${output}.f0.mat\n";
print MAT "OTU";
foreach $lib (sort keys %libhash){
    print MAT "\t$lib"
}
print MAT "\n";

foreach $seq (sort keys %indexhash){
    print MAT "$seq";
    foreach $lib (sort keys %libhash){
	if ($indexhash{$seq}{$lib}){
	    print MAT "\t$indexhash{$seq}{$lib}"; 
	} else {
	    print MAT "\t0";
	}
    }
    print MAT "\n";
}


