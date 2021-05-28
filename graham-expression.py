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
import os
import sys
import ProgramName
from BedReader import BedReader
from Interval import Interval
from Rex import Rex
rex=Rex()

# Note: Graham's peak coords are in hg38

#HIC="/data/reddylab/Tony/HiC/GC_Timecourse/merged_hics/iter3/0HR_merged.txt.hires.hic.fdr.0.01.loops.bedpe"
TADS="/data/reddylab/Tony/HiC/GC_Timecourse/merged_hics/iter3/TADs/0HR_merged.txt.hires.hic.tads.merged"
HIC=TADS
EXPRESSION_DATA="/data/reddylab/projects/GGR/results/rna_seq/differential_expression/iter0/edgeR/" # t05.vs.t00.protein_coding.txt
GENCODE="/data/reddylab/Reference_Data/Gencode/v22/gencode.v22.annotation.bed"
TSS="/home/bmajoros/PopSTARR/graham/tss.txt"

def loadTSS():
    byGene={}
    with open(TSS,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=4): continue
            (chr,pos,strand,gene)=fields
            pos=int(pos)
            rec=Interval(pos,pos+1)
            rec.chr=chr
            byGene[gene]=rec
    return byGene

def loadExpression(timepoint,tssHash):
    byChr={}
    filename=EXPRESSION_DATA+timepoint+".vs.t00.protein_coding.txt";
    IN=open(filename,"rt")
    IN.readline() # header
    for line in IN:
        fields=line.rstrip().split()
        (gene,logFC,logCPM,LR,P,Padj)=fields
        logFC=float(logFC)
        Padj=float(Padj)
        rec=tssHash.get(gene,None)
        if(rec is None): continue
        rec.logFC=logFC
        rec.Padj=Padj
        rec.type="gene"
        array=byChr.get(rec.chr,None)
        if(array is None): array=byChr[rec.chr]=[]
        array.append(rec)
    IN.close()
    return byChr

def loadHiC(expression):
    seen=set()
    with open(HIC,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            (chr1,begin1,end1,chr2,begin2,end2)=fields
            if(chr1!=chr2): continue
            array=expression.get(chr1,None)
            if(array is None): return
            anchor1=Interval(int(begin1),int(end1))
            anchor2=Interval(int(begin2),int(end2))
            anchor1.type=anchor2.type="anchor"
            anchor1.mate=anchor2
            anchor2.mate=anchor1
            anchor1.objects=[]; anchor2.objects=[]
            key1=chr1+" "+begin1+" "+end1
            key2=chr2+" "+begin2+" "+end2
            if(key1 not in seen): array.append(anchor1)
            if(key2 not in seen): array.append(anchor2)
            seen.add(key1); seen.add(key2)

def sortAndAssociate(expression):
    chroms=expression.keys()
    for chr in chroms:
        array=expression[chr]
        array.sort(key=lambda x: x.begin)
        #for elem in array: print(elem.begin,elem.end,elem.type,"\t")
        associate(array)

def associate(array):
    L=len(array)
    i=0
    while(i+1<L):
        thisRec=array[i]
        if(thisRec.type=="anchor"):
            scanLeft(thisRec,array,i-1)
            scanRight(thisRec,array,i+1)
        i+=1

def scanLeft(anchor,array,i):
    #print("scanLeft")
    while(i>=0):
        rec=array[i]
        if(rec.type=="anchor"): i-=1; continue
        if(not anchor.overlaps(rec)): return
        anchor.objects.append(rec)
        i-=1
        print(i)

def scanRight(anchor,array,i):
    #print("scanRight")
    L=len(array)
    while(i<L):
        rec=array[i]
        if(rec.type=="anchor"): i+=1; continue
        if(not anchor.overlaps(rec)): return
        anchor.objects.append(rec)
        i+=1

def addPeaks(bedDir,timepoint,direction,expression):
    records=BedReader.readAll(bedDir+"/"+timepoint+"_"+direction+".bed")
    for rec in records:
        array=expression.get(rec.chr,None)
        if(array is None): array=expression[rec.chr]=[]
        interval=rec.interval
        interval.begin=interval.intCenter()
        interval.end=interval.begin+1
        interval.type="peak"
        interval.dir=direction
        array.append(interval)

def loadPeaks(bedDir,timepoint,expression):
    addPeaks(bedDir,timepoint,"up",expression)
    addPeaks(bedDir,timepoint,"down",expression)
    addPeaks(bedDir,timepoint,"nonresp",expression)
    
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <bed-dir> <t8> \n")
(bedDir,timepoint)=sys.argv[1:]

print("loading TSSs...",file=sys.stderr,flush=True)
tssHash=loadTSS()

print("loading expression data...",file=sys.stderr,flush=True)
expression=loadExpression(timepoint,tssHash)

print("loading peaks...",file=sys.stderr,flush=True)
#upPeaks=BedReader.readAll(bedDir+"/"+timepoint+"_up.bed")
#downPeaks=BedReader.readAll(bedDir+"/"+timepoint+"_down.bed")
#samePeaks=BedReader.readAll(bedDir+"/"+timepoint+"_nonresp.bed")
loadPeaks(bedDir,timepoint,expression)

print("loading Hi-C data...",file=sys.stderr,flush=True)
loadHiC(expression)

# Now expression is a hash mapping chr to an array, each of which contains
# records of type "gene", "anchor", or "peak", and each peak has a 
# logFC and Padj

print("identifying associations...",file=sys.stderr,flush=True)
sortAndAssociate(expression)


