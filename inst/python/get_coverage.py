#!/usr/bin/python3

import sys, getopt
import pysam
import os
import numpy as np
import multiprocessing as mp
from datetime import datetime


def divide_chunks(lst, n):
    for i in range(0, len(lst), n): 
        yield lst[i:i + n]

def get_coverage(gps:None,bamfile:None,output:None,force:False,id:None):
   bamfile = pysam.AlignmentFile(bam, "rb")
   if gps!=None:
         for gpos in gps:
            chr=None
            start=None
            end=None
            ref="."
            alt="."
            gt="."
            gid="."
            af=0
            chr=gpos[0]
            if gpos[1]!=None:
               start=int(gpos[1])-1
               end=int(gpos[1])
            if gpos[2]!=None:
               ref=gpos[2]
            if gpos[3]!=None:
               alt=gpos[3]
            if gpos[4]!=None:
               gt=gpos[4]
            if gpos[5]!=None:
               gid=gpos[5]
   
            cov=bamfile.count_coverage(
               contig=chr,
               start=start,
               stop=end,
               quality_threshold=qt
            )

            baseCount={"A":sum(cov[0]),"C":sum(cov[1]),"G":sum(cov[2]),"T":sum(cov[3])}
            bTotal=baseCount["A"]+baseCount["C"]+baseCount["G"]+baseCount["T"]
            if not ref=="." and not alt==".":
               if (baseCount[ref]+baseCount[alt])>0:
                     af=baseCount[alt]/(baseCount[ref]+baseCount[alt])
               
            if output!=None:
               if not os.path.exists(output) or force:
                  f=open(output,"a")
                  if id!=None:
                     f.write("\t".join(map(str,[chr,end,ref,alt,gt,gid,baseCount["A"],baseCount["C"],baseCount["G"],baseCount["T"],bTotal,af,id]))+"\n")
                  else:
                     f.write("\t".join(map(str,[chr,end,ref,alt,gt,gid,baseCount["A"],baseCount["C"],baseCount["G"],baseCount["T"],bTotal,af]))+"\n")
                  f.close()
               else:
                  raise FileExistsError
            else:
               if id!=None:
                  print(chr,end,ref,alt,gt,gid,baseCount["A"],baseCount["C"],baseCount["G"],baseCount["T"],bTotal,af,id,sep="\t")
               else:
                  print(chr,end,ref,alt,gt,gid,baseCount["A"],baseCount["C"],baseCount["G"],baseCount["T"],bTotal,af,sep="\t")
         
if __name__ == "__main__":
   bam = None
   gpos=None
   region=[None]*2
   threads = mp.cpu_count() - 1
   qt= 15
   verbose=False
   output=None
   force=False
   id=None
   try:
      opts, args = getopt.getopt(sys.argv[1:],"h:b:g:c:p:q:t:v:o:f:i:",["bam=","gpos=","chr=","pos=","qt=","threads=","verbose=","output=","force=","id="])
   except getopt.GetoptError:
      print ('get_coverage.py -b <bam> -g [gpos] -c [chrom] -p [pos] -q [qt] -t [threads] -v [verbose] -o [output] -f [force] -i [id]')
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
         print ('get_coverage.py -b <bam> -g [gpos] -c [chrom] -p [pos] -q [qt] -t [threads] -v [verbose] -o [output] -f [force] -i [id]')
         sys.exit()
      elif opt in ("-b", "--bam"):
         bam = arg
      elif opt in ("-g", "--gpos"):
         gpos = arg
      elif opt in ("-c", "--chrom"):
         region[0] = arg
      elif opt in ("-p", "--pos"):
         region[1] = int(arg)
      elif opt in ("-t", "--threads"):
         if int(arg)<threads:
            threads = int(arg)
      elif opt in ("-q", "--qt"):
         qt = int(arg)
      elif opt in ("-v", "--verbose"):
         verbose = bool(arg)
      elif opt in ("-o", "--output"):
         output = arg
      elif opt in ("-f", "--force"):
         force=bool(arg)
      elif opt in ("-i", "--id"):
         id=arg

   gposcontent = []
   start=datetime.now()
   if verbose:
      print("• Start Time: ",start)
   
   if gpos!=None:
      if os.path.exists(gpos):
         with open(gpos) as f:
            for line in f:
                  gposcontent.append(list(line.strip().split()))
      else:
         raise ValueError("Path to genomic position file not found")
   elif region[0]!=None:
      gposcontent.append(region)
   else:
      raise ValueError("No genomic position was provided")

   chunk_size=int(len(gposcontent)/threads)
   chunks=list(divide_chunks(gposcontent,chunk_size))
   jobs=len(chunks)

   if output!=None:
      if not os.path.exists(output) or force:
         f=open(output,"w")
         if id!=None:
            f.write("chrom\tpos\tref\talt\tgt\tgid\tA\tC\tG\tT\tdepth\taf\tid\n")
         else:
            f.write("chrom\tpos\tref\talt\tgt\tgid\tA\tC\tG\tT\tdepth\taf\n")
         f.close()
      else:
         raise FileExistsError
   else:
      if id!=None:
         print("chrom\tpos\tref\talt\tgt\tgid\tA\tC\tG\tT\tdepth\taf\tid")
      else:
         print("chrom\tpos\tref\talt\tgt\tgid\tA\tC\tG\tT\tdepth\taf")


   if jobs==1 or threads==1:
      list(map(get_coverage,gposcontent,[bam],[output],[True],[id]))
   elif jobs>1:
      with mp.Pool(threads) as thread:
         list(thread.starmap(get_coverage,zip(chunks,[bam]*threads,[output]*threads,[True]*threads,[id]*threads)))
   end=datetime.now()
  
   if verbose:
      print("• End Time: ",end)
      print("• Run Time: ", end-start)
      print("• Jobs: ", jobs)
      print("• Threads: ", threads)

