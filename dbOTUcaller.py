#!/usr/bin/env python

''' Distribution-based clustering: version 2.0
Accurate method for creating operational taxonomic units (OTUs) from sequence data
Input requirements:
*input files in both OTU by library matrix and alignment files
*parameters such as the distance criteria, pvalue cutoff and abundance criteria that need to be satisfied in order to create an OTU
*output prefix name, such that it will be unique and log and err files can be created (currently only set to stdout)
'''
#needed to parse sequence files and alignments
#record version information here
version="DBC version 2.0 updated 11/18/14"

import Bio
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import numpy as np
import sys
import argparse
import csv
from datetime import datetime
#not sure if all of these things are needed
#somehow need to get chisquare 
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
base = importr('base')
chisqtest=rpy2.robjects.r("chisq.test")
import gc
import math

#Create custom defs
def JS(x,y): #Jensen-shannon divergence
   import warnings
   warnings.filterwarnings("ignore", category = RuntimeWarning)
   idx = np.where(x==0)
   x_pos = x.copy()
   x_pos[idx] = 1
   idy = np.where(y==0)
   y_pos = y.copy()
   y_pos[idy] = 1
   dx = x*np.log2(2*x_pos/(x_pos+y))
   dy = y*np.log2(2*y_pos/(x+y_pos))
   d = 0.5*sum(dx+dy)
   return d

def JSsqrt(x,y):
   d = JS(x,y)
   return d**0.5

def hamdist(str1, str2):
   """Count the # of differences between equal length strings str1 and str2"""
   #This doesn't exactly equal the clustalw values
   diffs = 0
   length = 0
   for ch1, ch2 in zip(str1, str2):
      if ch1 != ch2:
         diffs += 1
      if ch1.isalpha() and ch2.isalpha():
         length += 1
   #divide by the length
   percent=float(diffs)/float(length)
   return percent

def findlowest(dounaligned, str1, str2):
      (aligndist)=hamdist(str1,str2)
      if dounaligned:
         ustr1=str1.replace('-','')
         ustr2=str2.replace('-','')
         (ualigndist)=hamdist(ustr1, ustr2)
         return(min(aligndist,ualigndist))
      else:
         return(aligndist)

def runchisq(printverbose, i1, i2, OTUarray, reps):
   """Given the index of two ids (existing OTU i1 and candidate i2) and the OTU matrix without headers
   prepare them for chisquare analysis
   This means removing all values that are 0 in both
   running chi2_contingency, getting the residuals
   and testing if they are too small and need to be simulated"""
   #limit to only the two that are tested
   if printverbose: log.write("Starting runchisq with these idexes %d %s %d %s\n" % (i1,onlyOTUs[i1], i2, onlyOTUs[i2]))
   preparearray=OTUarray[ [i1,i2] ]
   #grab only the columns that have values in either
   testarray=preparearray[:, np.where((preparearray != 0).sum(axis=0) != 0)]
   L=[]
   for i in testarray:
       for j in i:
         for k in j:
            L.append(k)
   #THIS WORKS WITH L!!!!
   if printverbose:
      pstring=str(L)
      log.write(pstring)
      log.write("\n")

   if len(L):
      if len(L) == 2:
         return(1)
      else:
         m = robjects.r.matrix(robjects.IntVector(L), ncol=2)
         #figure out a way to clear the matrix in case there's an issue
         resnosim= chisqtest(m)
         all=0
         newall=float(all)
         g5=0
         newg5=float(g5)
         for i in range(0,len(resnosim[6])):
            newall+=1
            if i > 5:newg5+=1
            percent=newg5/newall
         if percent < 0.80:
            #check to see if the parent OTU i1 has a stored similated value
            if printverbose: log.write("percent too low (%f), simulate\n" % percent)
            res= chisqtest(m,simulate=True, B=reps)
            if printverbose: log.write("Done chisq.test with simulate %f\n" % res[2][0])
            gc.collect()
            return(res[2][0])
         else:
            if printverbose: log.write("Chisq test ok; return nonsim chisq: %f\n" % resnosim[2][0])
            gc.collect()
            return(resnosim[2][0])
      
      

def assignOTU(dounaligned, printverbose, distancecriteria, abundancecriteria, pvaluecutoff, JScorrect, existingOTUalignment, OTUtable1, onlyOTUs, alldata,  x, i):
   """Assign sequence into an existing set of OTUs
   
   distancecriteria, the max distance between sequences that can be merged into an OTU
   abundancecriteria, the fold-difference the existingOTU has to be over the tested OTU to be included
   pvalucutoff, the pvals lower than the cutoff will remain distinct OTUs
   JScorrect, a tuple of whether or not to correct JScorrect[0] and the value to use as the cutoff JScorrect[1]
   existingOTUalignment, an alignment dataset(?) with only the existing OTUs
   OTUtable1, with only ecological info to pass to pval
   onlyOTUs, with only the OTU names
   all_data, with the sum at the end
   x, the string from the alignment in rec = next((r for r in alignment if r.id == onlyOTUs[nindex]), None): str(rec.seq)
   i, the index of the sequence to be worked on which corresponds to OTUtable1, onlyOTUs and alldata"""

   #something about the existingOTUalignment, the OTUtale1 with only ecological info, onlyOTUs with the ids, 
   #the all_data with only the ecological info and the sum at the end
   #the x=str(alignment[0].seq) and the index of the OTU to search for

   #search the existing OTUs for the 10 closests
   if printverbose:log.write("This is the current index: %d %s %d\n" % (i, onlyOTUs[i], alldata[i][-1]))

   L=[]
   for y in existingOTUalignment:
      d=findlowest(dounaligned, x, str(y.seq))
      #append this to the current list of close sequences
      if d < distancecriteria:
         L.append((d, y.name))
      
   #sort the list by the first value of the tuple
   if L:
      sortedL=sorted(sorted(L, key=lambda tup: tup[1]))
      #now the first 10 values will be the ones you want
      #now test whether the distribution or abundance allows merge
      #for all of them change range(0,10) to range(0,len(sortedL))
      upper=len(sortedL)
      for j in range(0,upper):
         #this will be the id of the closest
         #this is the OTU
         OTUid=sortedL[j][1]
         #get the index of the OTUid
         jindex=np.where(onlyOTUs==OTUid)
         actualjindex=jindex[0][0]
         #does it fit the distance criteria
         if printverbose: log.write("This is the next most similar OTU rep: %s\n" % onlyOTUs[actualjindex])
         if sortedL[j][0] < distancecriteria:
            if printverbose: log.write("Passed distance criteria: %s\n" % onlyOTUs[actualjindex])
            #ok to continue
            #does it fit the distance criteria
            if alldata[actualjindex][-1] >= alldata[i][-1]*abundancecriteria:
               if printverbose: log.write("Passed abundance criteria: %s\n" % onlyOTUs[actualjindex])
               #it satisfies the abundance criteria
               #get the pvalue
               pval=runchisq(printverbose, i, actualjindex, OTUtable1, 10000)
               x=OTUtable1[i]
               y=OTUtable1[actualjindex]
               if pval < pvaluecutoff:
                  if printverbose: log.write("Did not pass pvaluecutoff\n")
                  if JScorrect[0]:
                     JSD=JS(x/sum(x),y/sum(y))
                     if printverbose: log.write("JSD %f %s %s\n" % (JSD, onlyOTUs[actualjindex], onlyOTUs[i]))
                     if JSD < JScorrect[1]:
                        if printverbose: log.write("Similar enough by JSD to merge even though did not pass pvalue cutoff: merge\n")
                        tup=(actualjindex, pval, sortedL[j][0])
                        return('merged', tup)
                  #if it's outside of the cutoff, its significant
                  #assign to another 
                  #print to log if verbose
                  continue
               else:
                  if printverbose: log.write("Passed pvalue: merge\n")
                  #this is the OTU to merge with
                  #print to the log if verbose
                  tup=(actualjindex, pval, sortedL[j][0])
                  return('merged', tup)
            else:
               #get the next closest one and print to the log if verbose
               #print to log if verbose
               if printverbose: log.write("Did not pass abundance criteria\n")
               continue
         else:
            if printverbose: log.write("Did not pass distance criteria\n")
            tup=('NA', 'NA', 'NA')
            return('not merged', tup)
      else:
         #this will be something to do if you don't find an OTU to merge into
         #not sure if this is necessary
         if printverbose: log.write("Exit without break\n")
         tup=('NA', 'NA', 'NA')
         return('not merged', tup)
   else:
      if printverbose: log.write("No existing OTU reps passed distance criteria\nExit without break\n")
      tup=('NA', 'NA', 'NA')
      return('not merged', tup)

def workthroughtable (logdict, dounaligned, printverbose, distancecriteria, abundancecriteria, pvaluecutoff, JScorrect, OTUtable1, OTUtable2, alignment, onlyOTUs):
   if printverbose: log.write("Start workingthroughtable\n")
   new_col = OTUtable1.sum(1)[...,None]
   all_data = np.append(OTUtable1, new_col, 1)

   #this sorts the whole array by the last column, but 
   rsortindex=all_data[:,-1].argsort()

   #now I can work through this from the reverse
   for nindex in reversed(rsortindex):
      if printverbose: log.write("Starting workthough on %s index %d\n" % (onlyOTUs[nindex], nindex))
      #work though the index values from most to least abundant
      #nindex is the index of the next most abundant
      if 'existingOTUalignment' in locals():
         if printverbose: log.write("OTUs exist, assignOTUs\n")
         #OTUs exist, see if it will fit into the existingOTU set
         #this tests both the genetic and ecological similarity
         if printverbose: log.write("nindex %d\n" % nindex)
         if printverbose: log.write("id %s\n" % onlyOTUs[nindex]) 
         #get the sequence record for the next OTU
         rec = next((r for r in alignment if r.id == onlyOTUs[nindex]), None)
         if printverbose: log.write("%s string\n" % str(rec.seq))
         
         if onlyOTUs[nindex] in logdict:
               if printverbose: log.write("Found %s in old log file\n" % onlyOTUs[nindex])
               if logdict[onlyOTUs[nindex]][0] == onlyOTUs[nindex]:
                  if printverbose: log.write("Was a parent in old log file %s\n" % onlyOTUs[nindex])
                  if printverbose: log.write("This is being created in existingOTUalignment %d %s\n" % (nindex, onlyOTUs[nindex]))
                  existingOTUalignment.append(rec)
                  #record the rep sequence in the list file dict
                  listdict[onlyOTUs[nindex]] = [onlyOTUs[nindex]]
                  #record the OTU table in the outtable dict
                  outtable[OTUtable2['OTU'][nindex]]=OTUtable1[nindex]
               else:
                  if printverbose: log.write("Was a child in old log file %s\n" % onlyOTUs[nindex])
                  if printverbose: log.write("this is being appended %d %s to %s\n" % (nindex, onlyOTUs[nindex], logdict[onlyOTUs[nindex]][0]))
                  print "This is the sequence %s and this is the parent %s\n" %(onlyOTUs[nindex], logdict[onlyOTUs[nindex]][0])
                  listdict[logdict[onlyOTUs[nindex]][0]].append(onlyOTUs[nindex])
                  outtable[logdict[onlyOTUs[nindex]][0]]=outtable[logdict[onlyOTUs[nindex]][0]] + OTUtable1[nindex]
         else:
            res=assignOTU(dounaligned, printverbose, distancecriteria, abundancecriteria, pvaluecutoff, JScorrect, existingOTUalignment, OTUtable1, onlyOTUs, all_data, str(rec.seq), nindex)
            if printverbose: log.write("Finished assignOTU\n")
            if res[0] == 'merged':
               if printverbose: log.write("The result is merged\n")
               #print out the the log that it's part of the mergeOTU
               pstring="Changefrom,%s,%s,Changeto,p.value,%f,Dist,%f,Done\n" % (onlyOTUs[nindex], onlyOTUs[res[1][0]], res[1][1], res[1][2])
               log.write(pstring)
               #record the OTU in the listdict
               if printverbose: log.write("this is being appended %d %s to %d %s\n" % (nindex, onlyOTUs[nindex], res[1][0], onlyOTUs[res[1][0]]))
               listdict[onlyOTUs[res[1][0]]].append(onlyOTUs[nindex])
               #add the new data to the outtable dict
               outtable[onlyOTUs[res[1][0]]]=outtable[onlyOTUs[res[1][0]]] + OTUtable1[nindex]
            else:
               #nothing came back, so create it as a new OTU
               if printverbose: log.write("the result is nothing\n")
               pstring="Parent,%s,Done\n" % (onlyOTUs[nindex])
               log.write(pstring)
               if printverbose: log.write("This is being created in existingOTUalignment %d %s\n" % (nindex, onlyOTUs[nindex]))
               rec = next((r for r in alignment if r.id == onlyOTUs[nindex]), None)
               existingOTUalignment.append(rec)
               #record the rep sequence in the list file dict
               listdict[onlyOTUs[nindex]] = [onlyOTUs[nindex]]
               #record the OTU table in the outtable dict
               outtable[OTUtable2['OTU'][nindex]]=OTUtable1[nindex]
      else:
         #theres nothing to merge with, it's a parent (it's the first one)
         if printverbose: log.write("No OTUs exist, begin\n")
         #check if it was a parent in the old log
         if onlyOTUs[nindex] in logdict:
            #This was found in the old log
            #Check to make sure it was a parent
            if logdict[onlyOTUs[nindex]][0] != onlyOTUs[nindex]:
               #Raise error
               raise AgreementError("Disagreement between current and prevous analysis\n")
         else:
            string="Parent,%s,Done\n" % (onlyOTUs[nindex])
            log.write(string)
         #create the existingOTUalignment to search through
         if printverbose: log.write("This is being created in existingOTUalignment %d %s\n" % (nindex, onlyOTUs[nindex]))
         existingOTUalignment= MultipleSeqAlignment([])
         rec = next((r for r in alignment if r.id == onlyOTUs[nindex]), None)
         existingOTUalignment.append(rec)
         #record the rep sequence in the list file dictionary
         listdict[onlyOTUs[nindex]] = [onlyOTUs[nindex]]
         #record the OTU counts in the outtable file dictionary
         outtable[OTUtable2['OTU'][nindex]]=OTUtable1[nindex]
         
def printresults(outlistfilename, outtablefilename, outfastafilename, listdict, outtable):
   log.write("Finished distribution-based clustering\n\n")
   
   #write out the OTU list
   log.write("WRITING RESULTS\nList file: %s\n" % outlistfilename)
   for x in listdict:
      outlisthand.write("\t".join(listdict[x]))
      outlisthand.write("\n")

   #write out the OTU table results
   log.write("OTU table: %s\n" % outtablefilename)
   headers=OTUtable2.dtype.names[:]
   outtablehand.write("\t".join(headers))
   outtablehand.write("\n")
   for x in outtable:
      outtablehand.write(x)
      outtablehand.write("\t")
      for y in outtable[x]:
         outtablehand.write("%f\t" % y)
      outtablehand.write("\n")

   #write out the rep fastas
   log.write("Fasta of OTU rep sequences: %s\n" % outfastafilename)
   for x in outtable:
      outfastahand.write(">%s\n" % x)
      rec = next((r for r in alignment if r.id == x), None)
      outfastahand.write("%s\n" % str(rec.seq))

   timestamp=str(datetime.now())
   string="\nEnding time: %s\n" % (timestamp)
   log.write(string)   

def almlabclue(printverbose):
   clue="WWSPPD"
   if printverbose: log.write("The next clue will remain an ENIGMA unless you know %s\n" % (clue))

def readoldlog(oldlog, printverbose):
   #This is where I'll put the code to read in the old log information
   logdict=dict()
   if printverbose: log.write("Using old log file in analysis: %s\n" %args.oldlog)
   with open(args.oldlog) as tsv:
      for line in csv.reader(tsv, delimiter=","):
         if line:
            #record old values and predetermined relationships
            if line[-1] == "Done":
               if line[0] == "Changefrom":
                  if printverbose: log.write("This was changed in old log file:child %s, parent %s \n" % (line[1], line[2]))
                  #do stuff to record this information 
                  #record the parent/child relationship in the listdict
                  logdict[line[1]]=(line[2], float(line[5]), float(line[7]))
               elif line[0] == "Parent":
                  if printverbose: log.write("This was a parent in old log file: parent %s\n" % line[1])
                  logdict[line[1]]=(line[1], 0, 0)
   return(logdict)

if __name__ == '__main__':
   parser = argparse.ArgumentParser(description='Create OTUs using ecological and genetic information (DBC version 2.0)')
   parser.add_argument('OTUtablefile', help='OTU table input')
   parser.add_argument('alignmentfile', help='alignment file input')
   parser.add_argument('output', help='unique prefix for output log, list, OTU table and fasta files')
   parser.add_argument('-d', '--dist_cutoff', type=float, default=0.1, help='maximum genetic variation allowed to be within the same population (i.e. OTU)')
   parser.add_argument('-k', '--k_fold', default=10, type=float, help='abundance criteria: existing OTU rep must have at least k-fold increase over the candidate sequence to be joined (default 10 for seq error only)')
   parser.add_argument('-p', '--pvalue', type=float, default=0.0001, help='pvalue cut-off: this could vary depending on the total number of libraries')
   parser.add_argument('-u', '--unaligned', action='store_true', help='use the unaligned sequence to correct alignment issues')
   parser.add_argument('-v', '--verbose', action='store_true', help='verbose option to work through method in log file')
   parser.add_argument('-s', '--split', type=str, help='input a list of the sequences clustered to the same percent as the distance cut-off to speed up the analysis')
   parser.add_argument('-j', '--useJS', type=float, help='Merge statistically significantly different sequences if below Jensen-Shannon cut-off?')
   parser.add_argument('-o', '--oldlog', type=str, help='Incorporate the results from an old log file?')
   args = parser.parse_args()

   #open output files to make writable
   #log file
   if args.oldlog:
      #I THINK I NEED TO ACTUALLY READ THE OLD FILE AFTER I LOAD ALL OF THE OTHER STUFF
      with open(args.oldlog) as tsv:
         for line in csv.reader(tsv, delimiter=","):
            if line:
               #record old values (you have to re-run with the same criteria and files to be valid)
               if line[0] == "Distance cutoff":
                  args.dist_cutoff=float(line[1])
               elif line[0] == "Abundance criteria":
                  args.k_fold=float(line[1])
               elif line[0] == "Pvalue cutoff":
                  args.pvalue=float(line[1])
               elif line[0] == "Input matfile":
                  args.OTUtablefile=line[1]
               elif line[0] == "Input alignmentfile":
                  args.alignmentfile=line[1]
               elif line[0] == "Output prefix":
                  args.output=line[1]
               elif line[0] == "Unaligned distance values":
                  if line[1]:
                     dounaligned=True
                  else:
                     dounaligned=False
                  args.unalign=line[1]
               elif line[0] == "verbose":
                  if int(line[1]):
                     printverbose=True
                  else:
                     printverbose=False
                  args.verbose=int(line[1])
               elif line[0] == "splitting analysis with":
                  args.split=line[1]
                  if args.split:
                     splitlist=[]
                     with open(args.split) as tsv:
                        for line in csv.reader(tsv, delimiter="\t"):
                           splitlist.append(line)
               elif line[0] == "Using JS":
                  args.useJS=line[1]

                  
      if args.useJS:
         JScorrect=(1, args.useJS)
      else:
         JScorrect=(0,)

      log =open(args.oldlog, 'a')
      log.write("START OVER\nReusing old log file in analysis,%s\n" %args.oldlog)

                  
   else:
      #open output files to make writable
      #log file
      logfilename="%s.log" % args.output
      log =open(logfilename, 'w')
      #list file                                                                                                                                                                                      
      #log some beginning information
      log.write("%s\n\nBeginning time: " % version)
      timestamp=str(datetime.now())
      string="%s\n" % (timestamp)
      log.write(string)
      log.write("\nSTARTING INPUT VALUES AND CONDITIONS\n")
      if args.unaligned:
         log.write("Unaligned distance values,1,will be evaluated\n")
         dounaligned=True
      else:
         log.write("Unaligned distance values,0,will not be evaluated\n")
         dounaligned=False
         
      if args.verbose: 
         log.write("verbose,1, is on\n")
         printverbose=True
      else:
         log.write("verbose,0, is off\n")
         printverbose=False
            
      if args.split:
         log.write("splitting analysis with,%s\n" % args.split)
         splitlist=[]
         with open(args.split) as tsv:
            for line in csv.reader(tsv, delimiter="\t"):
               splitlist.append(line)
   
      if args.useJS:
         log.write("Using JS,%f\n" % args.useJS)
         JScorrect=(1, args.useJS)
      else:
         JScorrect=(0,)

      log.write("Distance cutoff,%f\nAbundance criteria,%f\nPvalue cutoff,%f\n" % (args.dist_cutoff, args.k_fold, args.pvalue))
      log.write("Input matfile,%s\nInput alignmentfile,%s\nOutput prefix,%s\n" % (args.OTUtablefile, args.alignmentfile, args.output))


   #Open the files to gather information
   #list file                                                                                                                                                                                      
   outlistfilename="%s.list" % args.output
   outlisthand=open(outlistfilename, 'w')
   #fasta output file                                                                                                                                                                              
   outfastafilename="%s.fasta" % args.output
   outfastahand=open(outfastafilename, 'w')
   #mat file                                                                                                                                                                                       
   outtablefilename="%s.mat" % args.output
   outtablehand=open(outtablefilename, 'w')

   log.write("\nNOTES FROM METHOD\n")
   table = np.genfromtxt(args.OTUtablefile, comments="#")
   OTUtable1=table[1:,1:]
   OTUtable2 = np.genfromtxt(args.OTUtablefile, comments="#", names=True, dtype=None)
   alignment = AlignIO.read(args.alignmentfile, "fasta")
   onlyOTUs=OTUtable2['OTU']
   listdict=dict()
   outtable={}
   logdict=dict()
   if args.oldlog:
      logdict=readoldlog(args.oldlog, printverbose)

   if (args.split):
      if printverbose: log.write("Splitting into subclusters based on %s\n" %args.split)
      itno=0
      for cluster in splitlist:
         #clear the existing OTUs to make room for new ones
         itno += 1
         if printverbose: log.write("\n\nBeginning split no. %d\n" %itno)
         if 'existingOTUalignment' in locals():
            del existingOTUalignment
         listdict=dict()
         outtable={}
         #WITH OLD.LOG, SEE IF IT WAS DONE ALREADY
         #now do the workthrough for each cluster
         #begin to make subsets of OTUtables
         mask=np.ones(len(cluster), dtype=bool)
         templist=[]
         for isolate in cluster:
            #make a mask with only the index values for the isolates in this cluster with templist
            fakeindex=np.where(OTUtable2['OTU']==isolate)
            realindex=fakeindex[0][0]
            templist.append(realindex)

         #make the mask
         mask=[[templist]]
         #apply the mask to the tables to make subtables (although they return tuples and the first values are right
         subOTUtable1=OTUtable1[mask]
         subOTUtable2=OTUtable2[mask]
         subonlyOTUs=onlyOTUs[mask]
         workthroughtable (logdict, dounaligned, printverbose, args.dist_cutoff, args.k_fold, args.pvalue, JScorrect, subOTUtable1[0], subOTUtable2[0], alignment, subonlyOTUs[0])
         outlistfilename="%s.list" % args.output
         outtablefilename="%s.mat" % args.output
         outfastafilename="%s.fasta" % args.output
         printresults(outlistfilename, outtablefilename, outfastafilename, listdict, outtable)
   else:
      
      workthroughtable (logdict, dounaligned, printverbose, args.dist_cutoff, args.k_fold, args.pvalue, JScorrect, OTUtable1, OTUtable2, alignment, onlyOTUs)
      outlistfilename="%s.list" % args.output
      outtablefilename="%s.mat" % args.output
      outfastafilename="%s.fasta" % args.output
      printresults(outlistfilename, outtablefilename, outfastafilename, listdict, outtable)
      

      
