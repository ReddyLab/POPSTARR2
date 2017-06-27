#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
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
    def __init__(self,id,variants):
        self.id=id
        self.variants=variants

class Variant:
    def __init__(self,pos,id,allele):
        self.id=id
        self.pos=pos
        self.allele=allele

class Allele:
    def __init__(self,pos,seq):
        self.pos=pos
        self.seq=seq

def parseAttributes(fields,seq):
    variants=[]; alleles=[]
    for i in range(11,len(fields)):
        if(rex.find("Zs:Z:(\S+)",fields[i])):
            variants=parseVariants(rex[1],seq)
        #elif(rex.find("MD:Z:(\S+)",fields[i])):
        #    alleles=parseAlleles(rex[1])
    return (variants,alleles)

def parseVariants(text,seq):
    variants=[]
    fields=text.split(",")
    pos=0
    for field in fields:
        if(rex.find("(\d+)\|S\|(\S+)",field)):
            relpos=int(rex[1])
            pos+=relpos
            variant=Variant(pos,rex[2],seq[pos]) ###
            variants.append(variant)
            pos+=1
        elif(rex.find("(\d+)\|I\|(\S+)",field)):
            relpos=int(rex[1])
            pos+=relpos
        #    variant=Variant(pos,rex[2],seq[pos])
        #    variants.add(variant)
            pos+=1
        elif(rex.find("(\d+)\|D\|(\S+)",field)):
            relpos=int(rex[1])
            pos+=relpos
        #    variant=Variant(pos,rex[2],seq[pos])
        #    variants.add(variant)
            pos+=1
    return variants

#def parseAlleles(text):
#    alleles=[]
#    while(rex.find("^(\d+)(\D+)(\S*)",text)):
#        allele=Allele(int(rex[1]),rex[2])
#        alleles.append(allele)
#        text=rex[3]
#    return alleles

def getRead(IN):
    while(True):
        line=IN.readline()
        if(not line): return None
        if(rex.find("^@",line)): continue
        fields=line.rstrip().split()
        if(len(fields)<12): continue
        break
    (readID,flags,chrom,pos,x,cigar,equals,otherPos,d,seq,qual)=fields[:11]
    (variants,alleles)=parseAttributes(fields,seq)
    assignAllelesToVariants(variants,alleles)
    read=Read(readID,variants)
    return read

def assignAllelesToVariants(variants,alleles):
    for allele in alleles:
        for variant in variants:
            if(allele.pos==variant.pos):
                variant.allele=allele.seq
                break

def incAlleleCount(hash,id,allele):
    alleles=hash.get(id)
    if(alleles is None): alleles=hash[id]={}
    alleles[allele]=alleles.get(allele,0)+1

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.sam>\n")
(infile,)=sys.argv[1:]

# Process all reads
fragments=[]
lastRead=None
IN=open(infile,"rt")
while(True):
    read=getRead(IN)
    if(read is None): break
    if(lastRead and read.id==lastRead.id):
        for variant in read.variants:
            lastRead.variants.append(variant)
    elif(len(read.variants)>0): fragments.append(read)
    lastRead=read
IN.close()

# Count alleles
alleleCounts={}
for fragment in fragments:
    seen=set()
    for variant in fragment.variants:
        #key=variant.id+" "+variant.allele
        key=variant.id
        if(key in seen): continue
        incAlleleCount(alleleCounts,variant.id,variant.allele)
        seen.add(key)

# Generate output
variantsSorted=[]
for variantID in alleleCounts.keys(): variantsSorted.append(variantID)
variantsSorted.sort()
for variantID in variantsSorted:
    print(variantID,end="")
    alleles=alleleCounts[variantID]
    keysSorted=[]
    for allele in alleles.keys(): keysSorted.append(allele)
    keysSorted.sort()
    for allele in keysSorted:
        print("\t"+allele+"="+str(alleles[allele]),end="")
    print()

