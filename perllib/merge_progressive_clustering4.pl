#! /usr/bin/perl -w

die "Usage: lists,comma,del output\n" unless (@ARGV);

($lists, $output)=(@ARGV);
chomp ($lists);
chomp ($output);
die "Please follow comannd line options\n" unless ($output);

open (OUT, ">${output}") or die "Can't open ${output}\n";
(@files)=split (",", $lists);
foreach $file (@files){
    chomp ($file);
    open (IN, "<$file" ) or die "Can't open $file\n";
    while ($line =<IN>){
	chomp ($line);
	next unless ($line);
	($parent, @children)=split ("\t", $line);
	if (@children){
	    foreach $OTU (@children){
		$parenthash{$OTU}=$parent;
	    }

	}
	$allparenthash{$parent}++;
    }
}

#fix it so that only the upper level parents are represented as "parents"
foreach $child (keys %parenthash){
    $finished = ();
    $parentname = $parenthash{$child};
    until ($finished){
	if ($parenthash{$parentname} && $parenthash{$parentname} ne $parentname){
	    #means that the parent has a parent that is not itself
	    $parentname2 = $parenthash{$parentname};
	    $parentname=$parentname2;
	} else {
	    $finished++;
	}
    }
    $newparenthash{$child} = $parentname;
}

foreach $child (sort keys %newparenthash){
    $finalhash{$newparenthash{$child}}{$child}++;
    $finalallparents{$newparenthash{$child}}++;
}

foreach $parent (sort {$finalallparents{$b} <=> $finalallparents{$a}} keys %finalallparents){
    print  OUT "$parent";
    $printed{$parent}++;
    foreach $child (keys %{$finalhash{$parent}}){
	print OUT "\t$child";
	$printed{$child}++;
    }
    print OUT "\n";
}

foreach $parent (keys %allparenthash){
    print OUT "$parent\n" unless ($printed{$parent});
}
