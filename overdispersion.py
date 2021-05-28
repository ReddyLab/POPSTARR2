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
import os
import math
from Rex import Rex
rex=Rex()
from SummaryStats import SummaryStats
from Pipe import Pipe

VAR_SAMPLE_SIZE=1000
MIN_SAMPLE_SIZE=30
MIN_COUNT=1
NUM_QUANTILES=10
TRIM_PERCENTILE=0.999
MAX_RECORDS=10000000
INFILE="/home/bmajoros/PopSTARR/graham/overdispersion-fields.txt"
# Data originated from: /data/reddylab/gjohnson/whole_genome_STARRseq/wgss2_lipofectamine/hiseq/ase/iter_2/test_alleles/all

def load(filename,maxN):
    records=[]
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=4): exit()
            (DNAref,DNAalt,RNAref,RNAalt)=fields
            DNAref=int(DNAref); DNAalt=int(DNAalt)
            RNAref=int(RNAref); RNAalt=int(RNAalt)
            if(DNAref<MIN_COUNT or DNAalt<MIN_COUNT or
               RNAref<MIN_COUNT or RNAalt<MIN_COUNT): continue
            rec=[DNAref,DNAalt,RNAref,RNAalt]
            records.append(rec)
            if(len(records)>=maxN): break
    return records

def getQuantiles(records,index):
    X=[]
    for record in records: X.append(record[index])
    X.sort(key=lambda x: x)
    n=int(float(len(X))*TRIM_PERCENTILE)
    elemPerBin=float(n)/float(NUM_QUANTILES)
    if(elemPerBin<=0): raise Exception("elemPerBin="+str(elemPerBin))
    nextIndex=elemPerBin
    quantiles=[]
    while(nextIndex<n):
        x=X[int(nextIndex)]
        nq=len(quantiles)
        if(nq==0 or x>quantiles[nq-1]): quantiles.append(x)
        nextIndex+=elemPerBin
    quantiles.append(X[n-1])
    return quantiles

def inBounds(record,quantilesByVariate):
    for i in range(4):
        x=record[i]
        if(x<MIN_COUNT): return False
        quantiles=quantilesByVariate[i]
        n=len(quantiles)
        if(x>quantiles[n-1]): return False
    return True

def makeArray(dimensions):
    array=[]
    for i in range(dimensions[0]):
        row=[]
        for j in  range(dimensions[1]):
            col=[]
            for k in range(dimensions[2]):
                col.append([])
            row.append(col)
        array.append(row)
    return array

def lookupQuantile(x,quantiles):
    for i in range(len(quantiles)):
        if(x<=quantiles[i]): return i
    return len(quantiles)-1

def getArrayIndices(record,quantilesByVariate):
    point=[]
    for i in range(3):
        q=lookupQuantile(record[i],quantilesByVariate[i])
        point.append(q)
    return point

def processAll(records):
    # First, make a 3-dimensional array
    quantilesByVariate=[]
    dimensions=[]
    for index in range(4):
        quantiles=getQuantiles(records,index)
        quantilesByVariate.append(quantiles)
        dimensions.append(len(quantiles))
        #print(index,quantiles)
    cube=makeArray(dimensions)

    # Process the records, put counts into array
    for record in records:
        (DNAref,DNAalt,RNAref,RNAalt)=record
        if(not inBounds(record,quantilesByVariate)): continue
        point=getArrayIndices(record,quantilesByVariate)
        cube[point[0]][point[1]][point[2]].append(RNAalt)

    # Compute variance within each bin
    for i in range(dimensions[0]):
        for j in range(dimensions[1]):
            for k in range(dimensions[2]):
                samples=cube[i][j][k]
                n=len(samples)
                if(n<MIN_SAMPLE_SIZE): continue
                (mean,SD,min,max)=SummaryStats.roundedSummaryStats(samples)
                var=SD*SD
                (alpha,beta,N)=getMeanValues(i,j,k,quantilesByVariate)
                bbVar=betaBinomialVariance(alpha+1,beta+1,N)
                bbVar2=betaBinomialVariance2(alpha+1,beta+1,N)
                print(round(var,2),round(bbVar,2),round(bbVar2,2),
                      round(alpha,2),round(beta,2),round(N,2),i,j,k)

def geometricMean(a,b):
    sum=0
    for x in range(int(a),int(b)): sum+=math.log(x)
    #print(a,b,sum)
    return math.exp(sum/(b-a+1))

def getMeanValues(i,j,k,quantilesByVariate):
    point=[i,j,k]
    meanValues=[]
    for i in range(3):
        index=point[i]
        lower=1 if index==0 else quantilesByVariate[i][index-1]+1
        upper=quantilesByVariate[i][index]
        lower=float(lower); upper=float(upper)
        #mid=(lower+upper)/2.0
        mean=(upper*(upper+1)/2.0-lower*(lower-1)/2.0)/(upper-lower+1)
        #mean=geometricMean(lower,upper)
        meanValues.append(mean)
    return meanValues

def betaBinomialVariance(alpha,beta,n):
    var=n*alpha*beta*(alpha+beta+n)/((alpha+beta)*(alpha+beta)*(alpha+beta+1))
    return var

def betaBinomialVariance2(alpha,beta,n):
    pipe=Pipe("/home/bmajoros/src/util/beta-binomial -S "+str(VAR_SAMPLE_SIZE)+
              " "+str(n)+" "+str(alpha)+" "+str(alpha+beta))
    array=[]
    while(True):
        line=pipe.readline()
        if(line is None): break
        x=int(line.rstrip())
        array.append(x)
    (mean,SD,min,max)=SummaryStats.summaryStats(array)
    var=SD*SD
    if(var<0.01):
        print("VAR",array,"\n",alpha,beta,n,var)
        exit()
    return var

#=========================================================================
# main()
#=========================================================================
#if(len(sys.argv)!=2):
#    exit(ProgramName.get()+" <>\n")
#()=sys.argv[1:]

records=load(INFILE,MAX_RECORDS)
processAll(records)


