#! /usr/bin/perl -w
#
#	Use this program to make tab files for FileMaker
#

	die "Usage: seqs.fna OTUlist > Redirect output\n" unless (@ARGV);
	
	chomp (@ARGV);
	($seqsfile, $OTUfile) = (@ARGV);
chomp ($seqsfile);
chomp ($OTUfile);
$/ = ">";
open (IN, "<$seqsfile") or die "Can't open $seqsfile\n";
while ($line1 = <IN>){	
    chomp ($line1);
    next unless ($line1);
    ($header, @pieces) = split ("\n", $line1);
    ($unique, $derep)= split (" ", $header);
    $hash{$derep}{$unique}++;
}

close (IN);

$/="\n";
open (IN, "<${OTUfile}") or die "Can't open $OTUfile\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    (@OTUlist)=split ("\t", $line);
    ($OTUrep)=$OTUlist[0];
    print "$OTUrep";
    foreach $seq (@OTUlist){
	if ($hash{$seq}){
	    foreach $name (sort keys %{$hash{$seq}}){
		print "\t$name"
	    }
	} else {
	    die "Missing name seqs for $seq\n";
	}
    }
    print "\n";
}
close (IN);
