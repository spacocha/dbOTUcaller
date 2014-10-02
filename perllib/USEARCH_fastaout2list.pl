#! /usr/bin/perl -w
#
#	Use this program to make tab files for FileMaker
#

die "Usage: input_file output\n" unless (@ARGV);

chomp (@ARGV);
($file, $output) = (@ARGV);
$/ = ">";
open (OUT, ">${output}") or die "Can't open $file\n";
open (IN, "<$file") or die "Can't open $file\n";
while ($line1 = <IN>){	
    chomp ($line1);
    next unless ($line1);
    (@pieces) = split ("\n", $line1);
    ($info) = shift (@pieces);
    ($name)=$info=~/^(.+?);size=/;
    if ($info=~/up=otu/){
	#this sequence is a member of another OTU
	push (@{$uchash{$name}}, $name);
    } elsif ($info=~/up=member/){
	if ($info=~/top=(.+);size=.+;up=/){
	    ($parent)=$info=~/top=(.+?);size=.+;up=/;
	    push (@{$uchash{$parent}}, $name);
	} else {
	    die "Not sure where top is:$info\n";
	}
    } elsif ($info=~/up=chimera/){
	#it's a chimera, so don't bother with it
    } else {
	die "Not sure what to look for $line1\n";
    }		
}

close (IN);

foreach $seedname (sort keys %uchash){
    print OUT "$seedname";
    foreach $name (@{$uchash{$seedname}}){
	next if ($name eq $seedname);
	print OUT "\t$name";
    }
    print OUT "\n";
}

