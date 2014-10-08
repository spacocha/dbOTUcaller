#! /usr/bin/perl -w
#
#	Use this program to make tab files for FileMaker
#

	die "Use this to make a translation table between fasta files
Usage: fasta_files\n" unless (@ARGV);
	
chomp (@ARGV);

$/=">";
foreach $file (@ARGV){
    chomp ($file);
    open (IN, "<${file}") or die "Can't open $file\n";
    while ($line=<IN>){
	chomp ($line);
	next unless ($line);
	($name, @seqs)=split ("\n", $line);
	if ($name=~/^ID[0-9]+[A-z]\/ab=[0-9]+\//){
	    ($shrtname)=$name=~/^(ID[0-9]+[A-z])\/ab=[0-9]+\//;
	} else {
	    $shrtname=$name;
	}
	($seq) = join ("", @seqs);
	$hash{$file}{$seq}=$shrtname;
    }
    close (IN);
}
print "Name\tseq";
foreach $file (@ARGV){
    chomp ($file);
    print "\t$file";
}
print "\n";

foreach $file (@ARGV){
    chomp ($file);
    foreach $seq (sort keys %{$hash{$file}}){
	next if ($done{$seq});
	$done{$seq}++;
	$count++;
	print "MERGE${count}ID\t$seq";
	foreach $otherfile (@ARGV){
	    chomp ($otherfile);
	    if ($hash{$otherfile}{$seq}){
		print "\t$hash{$otherfile}{$seq}";
	    } else {
		print "\tNA";
	    }
	}
	print "\n";
    }
}
