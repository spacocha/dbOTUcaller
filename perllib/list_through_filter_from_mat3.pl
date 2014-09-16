#! /usr/bin/perl -w
#
#	Use this program as part of the eOTU pipeline
#

	die "Input the following:
full_list- file name of the complete OTU list
line_no- the line from the OTU list that will be evalulated
fasta- full fasta file with sequences and names that match those in the list
trans- trans file that relates the counts in the sequencing run with the unique fasta sequences
filter_number- the number of times the sequence must be seen in order to keep it
output_prefix- what the output files will begin with
comb_reps_file- file with two columns. The first column is the unique lib name, the second is the combreps lib name 
(multiple unique libs will be merged if they have the same combreps lib name)

Usage: full_list line_no fasta trans filter_number output_prefix comb_reps_file\n" unless (@ARGV);
	
	chomp (@ARGV);
	($full_list, $lineno,  $fullfasta, $fulltrans, $filterno, $output, $combreps) = (@ARGV);
die "list_through_filter_from_mat.pl: Please follow the command line arguments\n" unless ($output);

chomp ($full_list);
chomp ($lineno);
chomp ($fullfasta);
chomp ($fulltrans);
chomp ($filterno);
chomp ($output);
chomp ($combreps);
#read in the list file and capture the OTU on the right line according to the lineno
$curline=0;
open (IN, "<${full_list}") or die "Can't open ${full_list}\n";
while ($line1 = <IN>){
    chomp ($line1);
    $curline++;
    next unless ($curline == $lineno);
    (@pieces) = split ("\t", $line1);
    foreach $piece (@pieces){
	$otuhash{$piece}++;
    }
}
close (IN);

die "list_through_filter_from_mat.pl:Missing OTUhash\n" unless (%otuhash);
#reminder info: %otuhash contains all of the OTUs that we want to evaluate
#nothing else should be retained

#change the line breaks to > for fasta format
$/ = ">";

#read in the full fasta file
open (IN, "<$fullfasta") or die "Can't open $fullfasta\n";
while ($line1 = <IN>){	
    chomp ($line1);
    next unless ($line1);
    $sequence = ();
    (@pieces) = split ("\n", $line1);
    ($info) = shift (@pieces);
    #filter out all of the other entries
    next unless ($otuhash{$info});
    #merge the sequences from different lines
    ($sequence) = join ("", @pieces);
    #makes sure they don't contain the line breaks
    $sequence =~tr/\n//d;
    #store sequence infor for later processing
    $seqhash{$info}=$sequence;

}

close (IN);

die "Missing seqhash\n" unless (%seqhash);

#change back the line breaks
    $/="\n";
if ($combreps){
    open (IN, "<$combreps" ) or die "Can't open $combreps\n";
    while ($line =<IN>){
	chomp ($line);
	next unless ($line);
	($uniquelib, $comblib)=split ("\t", $line);
        $combhash{$uniquelib}=$comblib;
    }
    close (IN);

    die "Missing combhash after combreps\n" unless (%combhash);
}

#open in the full trans file
open (IN, "<$fulltrans" ) or die "Can't open $fulltrans\n";
while ($line =<IN>){
chomp ($line);
    next unless ($line);
    ($name, $OTU)=split ("\t", $line);
    if ($otuhash{$OTU}){
        if ($name=~/^.+_[0-9]+ /){
           ($lib)=$name=~/^(.+?)_[0-9]+ /;
        } elsif ($name=~/^.+_[0-9]+$/){
           ($lib)=$name=~/^(.+?)_[0-9]+$/;
        } else {
           die "Don't recognize the lib in $name\n";
        }
	if ($combreps){
	   if ($combhash{$lib}){
		$newlib=$combhash{$lib};
                $lib=$newlib;
	    } else {
		die "No value for $lib in combreps file $combreps: $combhash{$lib}\n";
	    }
	}
	$alllib{$lib}++;
	$OTU_lib_hash{$OTU}{$lib}++;
	$cutoffhash{$OTU}++;
    }
}
close (IN);

die "Missing OTU_lib_hash\n" unless (%OTU_lib_hash);
die "Missing cutoffhash\n" unless (%cutoffhash);

#make the mat file for use fa on the list and above the cutoff
open (MAT, ">${output}.mat") or die "Can't open ${output}.mat\n";
open (FA, ">${output}.fa") or die "Can't open ${output}.fa\n";
foreach $OTU (sort keys %OTU_lib_hash){
    print MAT "OTU";
    next unless ($cutoffhash{$OTU}>=$filterno);
    foreach $lib (sort keys %alllib){
	print MAT "\t$lib";
    }
    print MAT "\n";
    last;
}

foreach $OTU (sort keys %OTU_lib_hash){
    next unless ($cutoffhash{$OTU}>=$filterno);
    print MAT "$OTU";
    print FA ">$OTU\n$seqhash{$OTU}\n";
    foreach $lib (sort keys %alllib){
        if ($OTU_lib_hash{$OTU}{$lib}){
	    print MAT "\t$OTU_lib_hash{$OTU}{$lib}";
        } else {
	    print MAT "\t0";
        }
    }
    print MAT "\n";
}

close (MAT);
close (FA);

