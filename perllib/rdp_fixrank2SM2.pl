#! /usr/bin/perl -w

die "Usage: rdp_fixrank thresh mat\n" unless (@ARGV);
($rdp, $thresh, $mat) = (@ARGV);
chomp ($rdp);
chomp ($mat);
chomp ($thresh);
die "Please include command line arg\n" unless ($mat);

open (IN, "<$rdp" ) or die "Can't open $rdp\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    $line=~s/\"//g;
    (@pieces)=split ("\t", $line);
    if ($pieces[0]=~/^.+;size=[0-9]$/){
	($OTU)=$pieces[0]=~/^(.+);size=[0-9]$/;
    } elsif ($pieces[0]=~/^.+;count=[0-9]$/){
	($OTU)=$pieces[0]=~/^(.+);count=[0-9]$/;
    } elsif ($pieces[0]=~/^.+\/ab=[0-9]\/$/){
	($OTU)=$pieces[0]=~/^.+\/ab=[0-9]\/$/;
    } else {
	($OTU)=$pieces[0]=~/^(.+)$/;
    }
    $gothash{$OTU}++;
    if ($pieces[5]){
	$hash{$OTU}{"k"}=$pieces[5] if ($pieces[7]>$thresh);
    }
    if ($pieces[8]){
	$hash{$OTU}{"p"}=$pieces[8] if ($pieces[10]>$thresh);
    }
    if ($pieces[11]){
	$hash{$OTU}{"c"}=$pieces[11] if ($pieces[13]>$thresh);
    }
    if ($pieces[14]){
	$hash{$OTU}{"o"}=$pieces[14] if ($pieces[16]>$thresh);
    }
    if ($pieces[17]){
	$hash{$OTU}{"f"}=$pieces[17] if ($pieces[19]>$thresh);
    }
    if ($pieces[20]){
	$hash{$OTU}{"g"}=$pieces[20] if ($pieces[22]>$thresh);
    }
    if ($pieces[23]){
	$hash{$OTU}{"s"}=$pieces[23] if ($pieces[25]>$thresh);
    } 
}
close (IN);
$first=1;
@lineage=("k", "p", "c", "o", "f", "g", "s");
open (IN, "<$mat" ) or die "Can't open $mat\n";
while ($line =<IN>){
    chomp ($line);
    next unless ($line);
    (@pieces)=split ("\t", $line);
    if ($first){
	foreach $piece (@pieces){
	    print "$piece\t";
	}
	print "ConsensusLineage";
	$first=();
    } else {
	next unless ($gothash{$pieces[0]});
	foreach $piece (@pieces){
	    print "$piece\t";
	}
	foreach $rank (@lineage){
	    if ($hash{$pieces[0]}{$rank}){
		print "${rank}__$hash{$pieces[0]}{$rank}";
	    } else {
		print "${rank}__";
	    }
	    print ";" unless ($rank eq "s");
	}
    }
    print "\n";
}
close (IN);
