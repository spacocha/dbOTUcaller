#! /usr/bin/perl -w

die "Use this to remove - and . from fasta sequences
Usage: file output" unless (@ARGV);
($file, $output) = (@ARGV);
chomp ($file);

$/=">";
open (IN, "<$file" ) or die "Can't open $file\n";
open (OUT, ">${output}") or die "Can't open $file\n";
while ($line =<IN>){
    chomp ($line);
    ($name, @seqs)=split ("\n", $line);
    next unless ($name);
    ($sequence)=join ("", @seqs);
    $sequence =~ s/-//g;
    $sequence =~ s/\.//g;
    print OUT ">$name\n$sequence\n" if ($name && $sequence);
}

