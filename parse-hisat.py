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
from Rex import Rex
rex=Rex()

class Read:
    def __init__(self,id,variants,alleles):
        self.id=id
        self.variants=variants
        self.alleles=alleles

class Variant:
    def __init__(self,pos,id):
        self.id=id
        self.pos=pos

class Allele:
    def __init__(self,pos,seq):
        self.pos=pos
        self.seq=seq

def parseAttributes(fields):
    print(len(fields))
    variants=[]; alleles=[]
    for i in range(11,len(fields)):
        if(rex.find("Zs:Z:(\S+)",fields[i])):
            variants=parseVariants(rex[1])
        elif(rex.find("MD:Z:(\S+)",fields[i])):
            alleles=parseAlleles(rex[1])
    return (variants,alleles)

def parseVariants(text):
    variants=[]
    fields=text.split(",")
    for field in fields:
        if(rex.find("(\d+)\|S\|(\S+)",field)):
            variant=Variant(int(rex[1]),rex[2])
            variants.append(variant)
    return variants

def parseAlleles(text):
    alleles=[]
    while(rex.find("^(\d+)(\D+)",text)):
        allele=Allele(int(rex[1]),rex[2])
        alleles.append(allele)
    return alleles

def getRead(IN):
    while(True):
        line=IN.readline()
        if(not line): return None
        if(rex.find("^@",line)): continue
        fields=line.rstrip().split()
        if(len(fields)<12): continue
        break
    (readID,flags,chrom,pos,x,cigar,equals,otherPos,d,seq,qual)=fields[:11]
    (variants,alleles)=parseAttributes(fields)
    read=Read(readID,variants,alleles)
    return read

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.sam>\n")
(infile,)=sys.argv[1:]

IN=open(infile,"rt")
while(True):
    read=getRead(IN)
    if(read is None): break

    
IN.close()


