#! /usr/bin/perl -w
#
#


die "Usage: UC_file output\n" unless (@ARGV);
($file1, $output) = (@ARGV);

	die "Please follow command line args\n" unless ($output);
chomp ($file1);
chomp ($output);

open (OUT, ">${output}") or die "Can't open $output\n";
open (IN, "<$file1") or die "Can't open $file1\n";
$seedcount=0;
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    next if ($line=~/^\#/);
    ($name, $class, $id, $pars, $OTU)=split ("\t", $line);
    if ($name=~/\/ab=/){
	($newname)=$name=~/^(.+)\/ab=/;
	$name=$newname;
    } elsif ($name=~/;size=/){
	($newname)=$name=~/^(.+);size=/;
	$name=$newname;
    }
    if ($OTU=~/\/ab=/){
	($newOTUname)=$OTU=~/^(.+)\/ab=/;
	$OTU=$newOTUname;
    } elsif ($OTU=~/;size=/){
	($newname)=$OTU=~/^(.+);size=/;
	$OTU=$newname;
    }

    if ($class eq "otu"){
	push (@{$uchash{$name}}, $name);
    } elsif ($class eq "match"){
	push (@{$uchash{$OTU}}, $name);
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
