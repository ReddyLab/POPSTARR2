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

def load(filename):
    hash={}
    with open(filename,"rt") as IN:
        IN.readline() # header
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=12): continue
            (chrom,pos,variant,P,Padj,effect,DNAref,DNAalt,RNAref,
             RNAalt,ref,alt)=fields
            rnaFreq=float(RNAalt)/float(RNAref+RNAalt)
            dnaFreq=float(DNAalt)/float(DNAref+DNAalt)
            hash[variant]=(rnaFreq,dnaFreq)
    return hash

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <bowtie.txt> <hisat.txt>\n")
(bowtieFile,hisatFile)=sys.argv[1:]

bowtie=load(bowtieFile)
hisat=load(hisatFile)
keys=bowtie.keys()
rnaGreater=0; rnaLess=0; rnaSame=0
dnaGreater=0; dnaLess=0; dnaSame=0
for variant in keys:
    if(hisat.get(variant,None)==None): continue
    pair=bowtie[variant]
    (bowtieRNAfreq,bowtieDNAfreq)=pair
    pair=hisat[variant]
    (hisatRNAfreq,hisatDNAfreq)=pair
    if(hisatRNAfreq>bowtieRNAfreq): rnaGreater+=1
    elif(hisatRNAfreq<bowtieRNAfreq): rnaLess+=1
    elif(hisatRNAfreq==bowtieRNAfreq): rnaSame+=1
    if(hisatDNAfreq>bowtieDNAfreq): dnaGreater+=1
    elif(hisatDNAfreq<bowtieDNAfreq): dnaLess+=1
    elif(hisatDNAfreq==bowtieDNAfreq): dnaSame+=1

    if(hisatRNAfreq>bowtieRNAfreq): print("RNA\tINCREASE\t",variant)
    elif(hisatRNAfreq<bowtieRNAfreq): print("RNA\tDECREASE\t",variant)
    if(hisatDNAfreq>bowtieDNAfreq): print("DNA\tINCREASE\t",variant)
    elif(hisatDNAfreq<bowtieDNAfreq): print("DNA\tDECREASE\t",variant)

    #print(bowtieRNAfreq,hisatRNAfreq,sep="\t")
    #print(bowtieDNAfreq,hisatDNAfreq,sep="\t")
rnaTotal=float(rnaGreater+rnaSame+rnaLess)
dnaTotal=float(dnaGreater+dnaSame+dnaLess)
rnaGreater=float(rnaGreater)/rnaTotal
rnaLess=float(rnaLess)/rnaTotal
dnaGreater=float(dnaGreater)/dnaTotal
dnaLess=float(dnaLess)/dnaTotal
print("RNA: ",rnaGreater,"increase,",rnaSame,"same,",rnaLess,"decrease")
print("DNA: ",dnaGreater,"increase,",dnaSame,"same,",dnaLess,"decrease")




