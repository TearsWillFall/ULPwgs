#!/usr/bin/python

import sys, getopt
import pysam
import numpy as np

def main(argv):
   bam = ''
   bed= ''
   qt= 15
   try:
      opts, args = getopt.getopt(argv,"hi:b:q:",["ifile=","bfile=","qt="])
   except getopt.GetoptError:
      print ('get_coverage.py -i <bam> -b <bed> -q [qt]')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('get_coverage.py -i <bam> -b <bed> -q [qt]')
         sys.exit()
      elif opt in ("-i", "--ifile"):
         bam = arg
      elif opt in ("-b", "--bfile"):
         bed = arg
      elif opt in ("-q", "--qt"):
         qt = int(arg)

   bamfile = pysam.AlignmentFile(bam, "rb")
   bedcontent = []
   with open(bed)as f:
        for line in f:
            bedcontent.append(line.strip().split())
   for x in bedcontent:
       cov=bamfile.count_coverage(contig=x[0],start=int(x[2])-1,stop=int(x[2]),quality_threshold=qt)
       baseCount={"A":sum(cov[0]),"C":sum(cov[1]),"G":sum(cov[2]),"T":sum(cov[3])}
       bTotal=baseCount["A"]+baseCount["C"]+baseCount["G"]+baseCount["T"]
       refAllele=baseCount[x[4]]
       altAllele=baseCount[x[5]]
       try:
           vaf=(altAllele)/bTotal
       except ZeroDivisionError:
           vaf=0
       print (x[0],x[1],x[2],baseCount["A"],baseCount["C"],baseCount["G"],baseCount["T"],
       bTotal,x[4],refAllele,x[5],altAllele,vaf)

if __name__ == "__main__":
   main(sys.argv[1:])
