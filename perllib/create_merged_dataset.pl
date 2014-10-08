#! /usr/bin/perl -w

die "This is to merge QIIME formated matrix files picked from a reference. 
Library names and lineages can not have \"QIIME\" in them
merged_file: file where parents are denoted as \"was a parent in\"
translation table: translation that corresponds unqiue.fa across populations
matfiles_list: Input a file where each line is the name of a matfile to process
output_prefix: Prefix for the output files

Usage: merged_file translation_table matfiles_list output_prefix\n" unless (@ARGV);
($mergefile, $trans, $matfile, $output)=(@ARGV);
chomp ($matfile);
chomp ($mergefile);
chomp ($trans);
chomp ($output);

open (IN, "<${mergefile}") or die "Can't open $mergefile\n";
while ($line=<IN>){
    chomp ($line);
    next unless ($line=~/was a parent in/);
    ($parent)=$line=~/^(.+) was a parent in/;
    $parenthash{$parent}++;
}
close (IN);

@headers=();
open (IN, "<${trans}") or die "Can't open $trans\n";
open (OUT, ">${output}.merged.fa") or die "Can't opne ${output}.merged.fa\n";
while ($line=<IN>){
    chomp ($line);
    ($mergeid, $seq, @pieces)=split ("\t", $line);
    if (@headers){
	$i=0;
	$j=@pieces;
	until ($i >=$j){
	    ($sample)=$headers[$i]=~/^(.+)\/unique.fa/;
	    $translatehash{$sample}{$pieces[$i]}=$mergeid;
	    $i++;
	}
	print OUT ">${mergeid}\n$seq\n" if ($parenthash{$mergeid});

    } else {
	(@headers)=@pieces;
    }
}
close (IN);
close (OUT);

@headers=();
open (MAT, "<$matfile" ) or die "Can't open $matfile\n";
open (OUT, ">${output}.merged.mat") or die "Can't open ${output}.merged.mat\n";
while ($mat =<MAT>){
    chomp ($mat);
    next unless ($mat);
    @headers=();
    ($sample)=$mat=~/^(.+)\//;
    open (IN, "<$mat" ) or die "Can't open $mat\n";
    while ($line =<IN>){
	chomp ($line);
	next unless ($line);
	next if ($line=~/QIIME/);
	($OTU, @pieces)=split ("\t", $line);
	if (@headers){
	    if ($parenthash{$translatehash{$sample}{$OTU}}){
		$i=0;
		$j=@pieces;
		until ($i >=$j){
		    $value=$pieces[$i];
		    die "merge_mats.pl: $value (under $headers[$i]) is not a number. Use merge_mats_lins.pl when lineage is in the last field\n" unless ($value=~/^\d+/);
		    $header=$headers[$i];
		    $hash{$translatehash{$sample}{$OTU}}{$header}+=$value;
		    $allheaders{$header}++;
		    $i++;
		}
	    }
	} else {
	    (@headers)=@pieces;
	}
    }
    close (IN);
}
close (MAT);
foreach $OTU (sort keys %hash){
    print OUT "OTU";
    foreach $header (sort keys %allheaders){
	print OUT "\t$header";
    }
    #use the last name for the lineage
    print OUT "\t$linname\n";
    last;
}

foreach $OTU (sort keys %hash){
    print OUT "$OTU";
    foreach $header (sort keys %allheaders){
	if ($hash{$OTU}{$header}){
	    print OUT "\t$hash{$OTU}{$header}";
	} else {
	    print OUT "\t0";
	}
    }
    print OUT "\n";
}
