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
import ProgramName
import random
from scipy.stats import fisher_exact

def loadTest(filename):
    recs=[]
    with open(filename,"rt") as IN:
        #IN.readline() # header
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)<12): continue
            if(fields[0]=="chrom"): continue
            (chr,pos,variant,P,Padj,effect,DNAref,DNAalt,RNAref,RNAalt,ref,\
                 alt)=fields
            DNAref=int(DNAref); DNAalt=int(DNAalt)
            RNAref=int(RNAref); RNAalt=int(RNAalt)
            recs.append([variant,chr,pos,P,Padj,effect,DNAref,DNAalt,RNAref,RNAalt])
    return recs

def estimatePower(DNAref,DNAalt,rna,foldChange,iterations):
    # alt/(N-alt)=FC  =>  alt=FC*(N-alt)  =>  alt=FC*N-FC*alt  =>  
    #    =>  alt+FC*alt=FC*N  =>  (FC+1)*alt=FC*N  =>  alt=FC*N/(FC+1)
    altCount=foldChange*rna/(foldChange+1.0)
    refCount=rna-altCount
    altRate=altCount/float(rna) # rate at which alt allele produces transcripts
    refRate=refCount/float(rna) # rate at which ref allele produces transcripts
    dna=float(DNAref+DNAalt)
    dnaRefP=float(DNAref)/float(dna); dnaAltP=float(DNAalt)/float(dna)
    totalTests=0; positives=0
    for iter in range(iterations):
        sampledRef=0; sampledAlt=0
        while(sampledRef+sampledAlt<rna):
            if(random.uniform(0,1)<=dnaRefP):
                if(random.uniform(0,1)<=refRate): sampledRef+=1
            else:
                if(random.uniform(0,1)<=altRate): sampledAlt+=1
        (oddsRatio,P)=fisher_exact([[DNAref,DNAalt],[sampledRef,sampledAlt]])
        if(P<ALPHA): positives+=1
        totalTests+=1
    power=float(positives)/float(totalTests)
    return power


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+\
             " <test-alleles-output> <fold-change> <alpha> <iterations>\n")
(infile,foldChange,ALPHA,iterations)=sys.argv[1:]
foldChange=float(foldChange)
iterations=int(iterations)
ALPHA=float(ALPHA)

records=loadTest(infile)
for rec in records:
    (variant,chr,pos,P,Padj,effect,DNAref,DNAalt,RNAref,RNAalt)=rec
    effect=float(effect); P=float(P); Padj=float(Padj)
    rna=RNAref+RNAalt
    power=estimatePower(DNAref,DNAalt,rna,foldChange,iterations)
    print(variant,round(power,2),chr,pos,round(P,4),round(Padj,4),\
              round(effect,2),DNAref,DNAalt,RNAref,RNAalt,sep="\t")



