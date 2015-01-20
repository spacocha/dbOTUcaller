dbOTUcaller - new implementation of distribution-based clustering

Distribution-based clustering version 2.0

Introduction

This is a compact version of distribution-based clustering re-written in python and about 10-fold faster than the previous version with the same accuracy.

Please cite the original distribution-based clustering paper when using this program:
Appl Environ Microbiol. 2013 Nov;79(21):6593-603. doi: 10.1128/AEM.00342-13. Epub 2013 Aug 23.
Distribution-based clustering: using ecology to refine the operational taxonomic unit.
Preheim SP1, Perrotta AR, Martin-Platero AM, Gupta A, Alm EJ.

Outline:
1.) Requirements
2.) Quick-start guide to running dbOTUcaller
3.) Running dbOTUcaller from SmileTrain
4.) Additional tools
5.) Contact me

1.) Requirements:

These are the following verified programs and versions needed to run dbOTUcaller-
R (version 2.15.3)
python (version 2.7.3)
rpy2 (version 2.3.6)
numpy (version 1.5.1)
biopython

These are the required imports into python for the program. 
Test to make sure these are functioning in python:
import Bio
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np
import sys
import argparse
import csv
from datetime import datetime
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
base = importr('base')
chisqtest=rpy2.robjects.r("chisq.test")
import gc
import math

Errors in the import of any of these will cause the program to fail.


2.) Quick-start quide to running dbOTUcaller

Once all of the required programs and packages are installed and accessible by python, distribution-based clustering can be called with the following inputs
A.) OTU table input file
    This is a file with counts of 100% clusters (i.e. dereplicated sequences) in all of the libraries. The first row is headers with 'OTU' as the first column header, the first column is OTU ids.
B.) alignment file
    This is a file of all the dereplicated sequences aligned to a reference. Typically this is aligned to a dataset representing only the amplified portion of the 16S rRNA gene that is sequenced.
C.) output prefix
    A unique name that will be used as the prefix for the output files. Four output files will be generated, .log, .list, .fasta and .mat

There are many flags that can be used to alter the parameters of the clustering. Please type the following for a list of commands:
python dbOTUcaller.py --help

This generates the following output:
usage: dbOTUcaller.py [-h] [-d DIST_CUTOFF] [-k K_FOLD] [-p PVALUE] [-u] [-v]
                      [-s SPLIT] [-j USEJS]
                      OTUtablefile alignmentfile output

Create OTUs using ecological and genetic information (DBC version 2.0)

positional arguments:
  OTUtablefile          OTU table input
  alignmentfile         alignment file input
  output                unique prefix for output log, list, OTU table and
                        fasta files

optional arguments:
  -h, --help            show this help message and exit
  -d DIST_CUTOFF, --dist_cutoff DIST_CUTOFF
                        maximum genetic variation allowed to be within the
                        same population (i.e. OTU)
  -k K_FOLD, --k_fold K_FOLD
                        abundance criteria: existing OTU rep must have at
                        least k-fold increase over the candidate sequence to
                        be joined (default 10 for seq error only)
  -p PVALUE, --pvalue PVALUE
                        pvalue cut-off: this could vary depending on the total
                        number of libraries
  -u, --unaligned       use the unaligned sequence to correct alignment issues
  -v, --verbose         verbose option to work through method in log file
  -s SPLIT, --split SPLIT
                        input a list of the sequences clustered to the same
                        percent as the distance cut-off to speed up the
                        analysis
  -j USEJS, --useJS USEJS
                        Merge statistically significantly different sequences
                        if below Jensen-Shannon cut-off?


Example data can be found at the following website:
http://web.mit.edu/afs/athena.mit.edu/user/s/p/spacocha/Public/16S-workshop/16S-analysis-workshop-tar.tar.gz

Download and open the file above. Within the example_data2_res folder will be an example of the uniput OTU table file (unique.f0.good.mat), input alignment (unique.good.align) and output files with the previx unique.dbOTU (unique.dbOTU.log, unique.dbOTU.list, unique.dbOTU.fasta, unique.dbOTU.mat).


3.) Running dbOTUcaller from SmileTrain

Interfacing with SmileTrain, which is a python wrapper for USEARCH and custom perl and python scripts is a fast and easy way to go from raw data to processed files quickly and reliably.
SmileTrain can be accessed through GitHub:

https://github.com/almlab/SmileTrain

To install SmileTrain, follow the documentation on the wiki:

https://github.com/almlab/SmileTrain/wiki

Alternatively, a practical guide to running SmileTrain and dbOTUcaller can be found with the example datasets:
http://web.mit.edu/afs/athena.mit.edu/user/s/p/spacocha/Public/16S-workshop/16S-analysis-workshop-tar.tar.gz

In a file called Practical_notes_SmileTrain_dbOTU.docx

Follow the instructions within that document for quick installation guide and test data.

4.) Additional tools

A.) Merging different analysis

When running samples across different sequencing lanes, or combining separate analyses, you will have to use the following commands to merge datasets.
-This assumes you have run --k_fold 10 (to remove potential sequencing errors) on at least two different samples and you want to join them.
The following example assumes you have done two different dbOTU calling analyses in ../Reactor_and_Coupon/raw_data/dbOTU_dir and ../GroundWater/raw_data/dbOTUcaller_dir with the output files as indicated.

Merge the data together:
1.) Figure out how the sequences relate to each other- they must be the same length:
perl ~/lib/dbOTUcaller/perllib/translate_fasta_eOTU.pl ../Reactor_and_Coupon/raw_data/dbOTU_dir/unique.fa ../GroundWater/raw_data/dbOTUcaller_dir/unique.fa > translate_complete.txt

2.) Only keep sequences that are parents- the rest are considered as errors. Notice you use the unique.dbOTU.list from the output of the dbOTU analysis and the translate_complete.txt from step 1:
perl ~/lib/dbOTUcaller/perllib/merge_eOTU_results_parentsonly2.pl translate_complete.txt ../Reactor_and_Coupon/raw_data/dbOTU_dir/unique.dbOTU.list,../GroundWater/raw_data/dbOTUcaller_dir/unique.dbOTU.list > merged_output.txt

3.) make a text file which contains each mat file on a different line with ls:
ls ../Reactor_and_Coupon/raw_data/dbOTU_dir/unique.f0.mat ../GroundWater/raw_data/dbOTUcaller_dir/unique.f0.mat > matfiles

4.) Merge the different analyses into one using the output prefix "myoutput" and the output from steps 1, 2 and 3:
perl ~/lib/dbOTUcaller/perllib/create_merged_dataset.pl merged_output.txt translate_complete.txt  matfiles myoutput

files with myoutput.* can be used for additional analyses and should contain samples from both analyses. This can be done with two or more datasets

5.) Contact me:

Please create an issue to contact me about problems installing or running dbOTUcaller. I am happy to help, but documenting issues will help others who encounter the same issues.
