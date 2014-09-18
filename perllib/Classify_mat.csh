#! /bin/bash

#path to file where all of the data is
#load in the variables from the second command line input
BIN=~/bin/

#classify with RDP using the fixrank
java -Xmx600m -jar ${BIN}/rdp_classifier_2.3/rdp_classifier-2.3.jar -q ${1} -o ${1}.rdp -f fixrank
#add classifications onto mat
perl ${BIN}/rdp_fixrank2SM2.pl ${1}.rdp 0.5 ${2} > ${2}.lin