#! /usr/bin/perl -w

die "Use this program to evaulate the find_seq_errors.pl 
Be in a folder where you ran make_it_go

Usage: vars_files output_prefix\n" unless (@ARGV);
($vars, $output)=(@ARGV);
chomp ($vars);
chomp ($output);
die "Please follow command line options\n" unless ($output);
$verbose=0;

print LOG "Creating output\n" if ($verobse);
open (LOG, ">${output}.log") or die "Can't open ${output}.log\n";

print LOG "Reading in variables\n" if ($verbose);
open (IN, "<$vars" ) or die "Can't open $vars\n";
while ($line=<IN>){
    chomp ($line);
    #get the PCLUSTOUTPUT name so you can check to see how many there are supposed to be

    if ($line=~/UNIQUE=/){
	($UNIQUE)=$line=~/UNIQUE=(.+)$/;
	print LOG "UNIQUE $UNIQUE\n";
    }
    if ($line=~/FNUMBER=/){
	($FNUMBER)=$line=~/FNUMBER=(.+)$/;
	print LOG "FNUMBER $FNUMBER\n";
    }
}

close (IN);
($PCLUSTOUTPUT)="${UNIQUE}.PC.final.list";

print LOG "Close variable file\n" if ($verbose);
#look for the number of lines in the PCLUSTOUTPUT


if ("-e $PCLUSTOUTPUT"){
    print LOG "$PCLUSTOUTPUT exists. Looking for lines\n" if ($verbose);
    $totlines=`cat $PCLUSTOUTPUT | wc -l`;
    chomp ($totlines);
    print LOG "Looking for $totlines of each file type\n";
} else {
    die "EVAL WARNING: Did not generate the PCLUSTOUTPUT $PCLUSTOUTPUT
EVAL ACTION: Start from the ProgressiveClustering program to generate this file\n";
}
#Looking for the error files
print LOG "Looking for the error files\n" if ($verbose);
if ("-e ${UNIQUE}*.f${FNUMBER}.err"){
    print LOG "Starting to look for error files\n" if ($verbose);
    @errfiles = `ls ${UNIQUE}*.f${FNUMBER}.err`;
    chomp (@errfiles);
    $numbererrs = @errfiles;
    chomp ($numbererrs);
    print LOG "Found $numbererrs error files\n" if ($verbose);
   if ($numbererrs == $totlines){
        print LOG "All of the $totlines of error files (${UNIQUE}*.f${FNUMBER}.err) were created\n";
    } else {
        print LOG "EVAL WARNING: Missing some error files: $numbererrs out of $totlines were created
EVAL ACTION: Rerun eOTU_core3.csh with the following lines:\n";
	$filetype = "error files";
	print LOG "Beginning prono search\n" if ($verbose);
	($result) = printmissingno ($totlines, $filetype);
    }
    print LOG "Evaluating how many of the $numbererrs are finished\n";
    foreach $f (@errfiles){
	chomp ($f);
	print LOG "Working on f $f\n" if ($verbose);
	$finished= `cat $f | grep "Finished"`;
	print LOG "Cat complete $f\n" if ($verbose); 
	chomp ($finished);
	$unfin=0;
	unless ($finished){
	    print LOG "EVAL WARNING:${f} did not finish\n";
# Finish this so that it will resubmit the error programs
	    ($tnum)=$f=~/\.([0-9])\.process\.f${FNUMBER}\.err/;
	    die "Missing t num $tnum $f ${FNUMBER} print \n" unless ($tnum);
	    $unfin++;
	}
    }
    print LOG "All of the error files finished\n" unless ($unfin);
} else {
    print LOG "EVAL WARNING: ${UNIQUE}*.f${FNUMBER}.err files do not exist
EVAL ACTION: Rerun eOTU_core program\n";
}

#look too see whether files were generated
if ("-e ${UNIQUE}*.process.f{$FNUMBER}.log"){
    @logfiles = `ls ${UNIQUE}*.process.f${FNUMBER}.log`;
    $numberlogs = @logfiles;
    if ($numberlogs == $totlines){
	print LOG "All of the $totlines of log files (${UNIQUE}*.process.f${FNUMBER}.log) were created\n";
    } else {
	print LOG "EVAL WARNING: Missing some log files: $numberlogs out of $totlines were created
EVAL ACTION: Rerun eOTU_core3.csh with the following lines\n";
	$filetype="log";
	($result) = printmissingno (@logfiles, $totlines, $filetype);
    }
    $warning=();
    foreach $file (@logfiles){
	$result=`cat ${file}`;
	if ($result=~/WARN/){
	    print LOG "There was at least one  WARNING in ${file}, check to see what the warning is\n";
	    $warning++;
	}
    }
    print LOG "There were no WARNINGs\n" unless ($warning);
} else {
    print LOG "EVAL WARNING: ${UNIQUE}*.process.f${FNUMBER}.log files do not exist
EVAL ACTION: Rerun the eOTU_core program\n";
    
}

sub printmissingno {

    my ($total, $filetype)= @_;
    my ($f);
    my ($i);
    my (%gothash);
    my ($done);
    %gothash=();
    $done=0;

    if ($filetype eq "error files"){
	foreach $f (@errfiles){
	    chomp ($f);
	    next unless ($f=~/${UNIQUE}.*.process/);
	    ($process)=${f}=~/${UNIQUE}.([0-9]+).process/;
	    die "F: $f Process: $process\n" unless ($process);
	    $gothash{$process}++;
	}
	$i=1;
	if ($total){
	    until ($i >=$total){
		print LOG "$i\t$filetype\n" unless ($gothash{$i});
		$i++;
	    }
	    $done++;
	    return ($done);
	} else {
	    die "Missing total $total\n";
	}
    }
}
