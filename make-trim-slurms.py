#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2017 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import os
import ProgramName
from SlurmWriter import SlurmWriter
from Rex import Rex
rex=Rex()

ROOT="/home/bmajoros/PopSTARR/graham"
MEM=50000
NICE=500
jobName="TRIM"
maxParallel=1000
THREADS=31
TRIMMOMATIC="java -jar /data/reddylab/software/Trimmomatic-0.33/Trimmomatic-0.33/trimmomatic-0.33.jar PE"

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <adapters.fasta> <fastq-in> <fastq-out> <full-path-to-slurms>\n")
(adaptersFasta,fastqIn,fastqOut,slurmDir)=sys.argv[1:]

files=os.listdir(fastqIn)
writer=SlurmWriter()
for file in files:
    if(not rex.find("(.*[_-])R1([_-].*)\.fastq.gz",file)): continue
    file1=file
    file2=rex[1]+"R2"+rex[2]+".fastq.gz"
    cmd=TRIMMOMATIC+" -threads "+str(THREADS)+" -phred33 "+\
        fastqIn+"/"+file1+" "+fastqIn+"/"+file2+" "+\
        fastqOut+"/"+rex[1]+"_FWD_paired.fq.gz "+\
        fastqOut+"/"+rex[1]+"_FWD_unpaired.fq.gz "+\
        fastqOut+"/"+rex[1]+"_REV_paired.fq.gz "+\
        fastqOut+"/"+rex[1]+"_REV_unpaired.fq.gz "+\
        "ILLUMINACLIP:"+adaptersFasta+\
        ":2:30:10:8:TRUE HEADCROP:1 LEADING:30 TRAILING:30 "+\
        "SLIDINGWINDOW:4:15 MINLEN:36"
    writer.addCommand("cd "+ROOT+"\n"+cmd)
writer.nice(NICE) # turns on "nice" (sets it to 100 by default)
writer.mem(MEM)
writer.threads(THREADS)
writer.setQueue("new,all")
writer.writeArrayScript(slurmDir,jobName,maxParallel,
                        "#SBATCH --exclude=x2-01-1,x2-01-2,x2-01-3,x2-01-4,x2-02-1,x2-02-2,x2-02-3,x2-02-4,x2-03-1 ")




