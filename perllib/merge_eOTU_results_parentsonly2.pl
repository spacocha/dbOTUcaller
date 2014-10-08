#! /usr/bin/perl
#
#

	die "Use this program to merge the results of the distribution-based clustering algorthim run on subset of total data
Make the  translation table using translate_fasta_eOTU.pl
Make the err.list files using err2list2.pl
Order the files in the same order as the translation table
Usage: translation_table comma,sep,err,files > redirect output\n" unless (@ARGV);
	($trans, $errs) = (@ARGV);

	die "Please follow command line args\n" unless ($errs);
chomp ($trans);
chomp ($errs);

(@errfiles)=split (",", $errs);

open (IN, "<$trans") or die "Can't open $trans\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line);
    ($mergeid, $seq, @pieces)=split ("\t", $line);

    if (@headers){
	$j = @headers;
	$i =0;
	until ($i >=$j){
	    $transhash{$errfiles[$i]}{$pieces[$i]}=$mergeid;
	    $i++;
	}
    } else {
	(@headers)=@pieces;
	$j = @headers;
	$i =0;
	until ($i >=$j){
	    print "$errfiles[$i]\t$headers[$i]\n";
            $i++;
	}
    }

}
close (IN);

foreach $file (@errfiles){
    chomp ($file);
    open (IN, "<${file}") or die "Can't open $file\n";
    while ($line=<IN>){
	chomp ($line);
	next unless ($line);
	($parent, @pieces)=split ("\t", $line);
	if ($transhash{$file}{$parent}){
	    #record as a parent in this file
	    $parenthash{$transhash{$file}{$parent}}{$file}++;
	    foreach $child (@pieces){
		if ($transhash{$file}{$child}){
		    $childhash{$transhash{$file}{$child}}{$file}=$transhash{$file}{$parent};
		} else {
		    die "Missing translation for $child\n";
		}
	    }
	} else {
	    die "Missing translation for $parent\n";
	}
    }
}

foreach $parent (sort keys %parenthash){
    print "$parent";
    #it was obviously a parent, since it's in the parent hash
    print " was a parent in:";
    foreach $file (sort keys %{$parenthash{$parent}}){
	print "$file;";
    }
    print "\n";
}


sub revcomp{
    (@pieces) = split ("", $rc_seq);
    #make reverse complement
    $j = @pieces;
    $j-=1;
    $seq = ();
    until ($j < 0){
	if ($pieces[$j]){
	    if ($pieces[$j] eq "A"){
		$seq = "$seq"."T";
	    } elsif ($pieces[$j] eq "T"){
		$seq = "$seq"."A";
	    } elsif ($pieces[$j] eq "C"){
		$seq = "$seq"."G";
	    } elsif ($pieces[$j] eq "G"){
		$seq = "$seq"."C";
	    } elsif ($pieces[$j] eq "N"){
		$seq = "$seq"."N";
	    } else {
		die "$pieces[$j]\n";
	    }
	} else {
	    die "NO j $j\n";
	}
	$j--;
    }
    return ($seq);
}
