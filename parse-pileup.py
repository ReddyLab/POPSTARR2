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
import gzip
from Rex import Rex
rex=Rex()

ACGT=set(["A","C","G","T"])
QUAL="!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ"
QUALITY=[-1]*256
for i in range(len(QUAL)): QUALITY[ord(QUAL[i])]=i

class Variant:
    def __init__(self,pos,id,ref,alt):
        self.id=id
        self.pos=pos
        self.ref=ref
        self.alt=alt

def removeIndels(bases):
    while(rex.find("(\S*)[+-](\d+)(\S*)",bases)):
        left=rex[1]; indelLen=int(rex[2]); right=rex[3]
        bases=left+right[indelLen:]
    while(rex.find("(\S*)\^.(\S*)",bases)):
        bases=rex[1]+rex[2]
    while(rex.find("(\S*)\$(\S*)",bases)):
        bases=rex[1]+rex[2]
    return bases

def parseBases(bases,qual):
    bases=removeIndels(bases)
    bases=bases.upper()
    counts={}
    #for c in bases:
    #    if(c in ACGT): counts[c]=counts.get(c,0)+1
    index=0
    for i in range(len(bases)):
        c=bases[i]
        if(c in ACGT):
            if(index>=len(qual)): exit(bases+"\n"+qual+"\n"+str(index)+"\t"+c)
            q=qual[index]
            #print(c,q,ord(q),sep="\t")
            if(QUALITY[ord(q)]>=MIN_QUAL): counts[c]=counts.get(c,0)+1
            index+=1
    return counts

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=4):
    exit(ProgramName.get()+" <sorted.vcf.gz> <pileup.txt> <min-quality>\n")
(vcf,pileup,MIN_QUAL)=sys.argv[1:]
MIN_QUAL=int(MIN_QUAL)

# Process the VCF file to get the ref and alt alleles
variants={}
for line in gzip.open(vcf):
    line=line.decode("utf-8")
    if(len(line)>0 and line[0]=="#"): continue
    fields=line.rstrip().split()
    if(len(fields)<9): continue
    (chr,pos,id,ref,alt,x,Pass,flags,GT)=fields[:9]
    altFields=alt.split(",")
    alt=altFields[0] # Keeping only the first alternate allele
    #if(id=="."): continue
    if(id=="."): id=chr+"@"+pos
    if(variants.get(chr,None) is None): variants[chr]=[]
    variants[chr].append(Variant(int(pos),id,ref,alt))
    #variants[chr].append(Variant(int(pos)-1,id,ref,alt)) # 8/21 added -1

# Process the pileup file
prevChr=None
variantsOnChr=None
nextVariant=None
for line in open(pileup,"rt"):
    fields=line.split()
    if(len(fields)<6): continue
    (chr,pos,N,total,seq,qual)=fields
    pos=int(pos)
    if(prevChr is None or chr!=prevChr):
        variantsOnChr=variants.get(chr,None)
        if(variantsOnChr is None): continue
        prevChr=chr
        nextVariant=0
    while(nextVariant<len(variantsOnChr) and
          variantsOnChr[nextVariant].pos<pos-1):
        nextVariant+=1
    if(nextVariant>=len(variantsOnChr)): continue
    variant=variantsOnChr[nextVariant]
    if(variant.pos!=pos): continue
    if(len(variant.ref)>1 or len(variant.alt)>1): continue # indel
    counts=parseBases(seq,qual)
    #print(line)
    refCount=counts.get(variant.ref,0)
    altCount=counts.get(variant.alt,0)
    print(variant.id,chr,variant.pos,variant.ref,variant.alt,refCount,altCount,sep="\t")


