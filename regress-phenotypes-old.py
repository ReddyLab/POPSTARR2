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
import math
import numpy as np
from sklearn import linear_model
import statsmodels.api as sm
from scipy import stats
from Rex import Rex
rex=Rex()

NUM_COVARIATES=11
SKIP_NA=True
BASE="/home/bmajoros/PopSTARR/sarah"
#EFFECTS=BASE+"/test-lucif-aug28-chr10.txt"
EFFECTS=BASE+"/test-sarah-beta-2sided.txt"
VCF=BASE+"/vcf"
PHENOTYPES=BASE+"/phenotypes.txt"
CHROMS=("chr10","chr11","chr21","chr5","chr8")

CENTERS=("A","E","F","K","N","P","Q")
CENTER_CODES={}
nextCode=0
for center in CENTERS:
    CENTER_CODES[center]=nextCode
    nextCode+=1

def integrateNAs(array):
    if(len(array)==0): return
    L=len(array[0])
    for i in range(5,L):
        sum=0; n=0
        for rec in array:
            x=rec[i]
            if(x!="NA" and rex.find("\d",x)):
                sum+=float(x)
                n+=1.0
        mean=sum/n
        for rec in array:
            if(rec[i]=="NA"): rec[i]=mean
        

def loadPhenotypes(filename):
    hash={}
    with open(filename,"rt") as IN:
        IN.readline()
        for line in IN:
            if(SKIP_NA and rex.find("NA",line)): continue
            fields=line.rstrip().split()
            (ID,otherID,FPG,PG1hr,PG2hr,FCP,CP1hr,pc1,pc2,age,gestationalAge,
             centerID,children,BMI,height,BP,smoker,drinker,ethnicity,
             ancestry)=fields
            FPG=math.log(float(FPG),10.0)
            FCP=math.log(float(FCP),10.0)
            PG1hr=math.sqrt(float(PG1hr))
            PG2hr=math.sqrt(float(PG2hr))
            CP1hr=math.sqrt(float(CP1hr))
            age=float(age); gestationalAge=float(gestationalAge)
            children=int(children); BMI=float(BMI); height=float(height)
            BP=float(BP); smoker=int(smoker); drinker=int(drinker)
            ethnicity=int(ethnicity)
            centerID=CENTER_CODES[centerID]
            rec=[FPG,PG1hr,PG2hr,FCP,CP1hr,age,gestationalAge,centerID,
                 children,BMI,height,BP,smoker,drinker,ethnicity]
            hash[ID]=rec
    return hash
#200848338@1097100918    P0376   70.19999695     117     84.59999847     0.5     6.6     0.006554579     0.011527396     25.1    28.1    P       1       23.73866272     168     72.33333588     0       0       1       EU

def loadEffects(filename):
    hash={}
    with open(filename,"rt") as IN:
        IN.readline()
        for line in IN:
            fields=line.rstrip().split()
            #(chrom,pos,variant,P,Padj,effect,DNAref,DNAalt,RNAref,RNAalt,
             #refAllele,altAllele)=fields
            (chrom,pos,variant,power,P,Padj,effect,DNAref,DNAalt,RNAref,RNAalt,
             refAllele,altAllele)=fields
            effect=float(effect)
            hash[variant]=[effect,refAllele,altAllele,P,Padj]
    return hash
#chr10	70827089	rs4745974	0.21501551605	1.0	0.31382978723404253	2	94	4	59	A      G

def processVariant(variant,VCFref,VCFalt,IDs,genotypes,phenotypes,effects,
                   phenotypeID):
    effect=effects.get(variant,None)
    if(effect is None): return
    (effectSize,refAllele,altAllele,P,Padj)=effect
    if(refAllele==VCFalt):
        if(effectSize==0.0): return
        effectSize=1.0/effectSize
    numIDs=len(IDs)
    points=[]
    for i in range(numIDs):
        ID=IDs[i]
        phenotypeRec=phenotypes.get(ID,None)
        if(phenotypeRec is None): continue
        phenotype=phenotypeRec[phenotypeID]
        genotype=genotypes[i]
        if(not rex.find("(\d+)\|(\d+)",genotype)): return
        allele1=int(rex[1]); allele2=int(rex[2])
        numAltAlleles=allele1+allele2
        totalEffect=numAltAlleles*effectSize
        point=[phenotype,totalEffect]
        point.extend(phenotypeRec[5:])
        points.append(point)
    regress(points,variant)

def regress(points,variantID):
    X=[]; Y=[]
    for point in points:
        #X.extend(point[1:])
        X.append(point[1])
        Y.append(point[0])
    X=np.array(X)
    Y=np.array(Y)
    #X=np.reshape(X,(-1,NUM_COVARIATES))
    X=np.reshape(X,(-1,1))
    Y=np.reshape(Y,(-1,1))
    regr = linear_model.LinearRegression()
    regr.fit(X,Y)
    #print('Coefficients: \n', regr.coef_)
    #print("Mean squared error: %.2f" % np.mean((regr.predict(X)-Y)**2))
    #print('r-squared: %.2f' % regr.score(X,Y))
    X2 = sm.add_constant(X)
    est = sm.OLS(Y, X2)
    est2 = est.fit()
    #print(est2.summary())
    #print(dir(est2))
    #print("P=",est2.f_pvalue)
    P=est2.f_pvalue
    if(P<0.05):
        MSE=np.mean((regr.predict(X)-Y)**2)
        r2=regr.score(X,Y)
        if(r2<0.04): return
        print(P,variantID,r2,sep="\t",flush=True)


#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <phenotype:1-5>\n")
(phenotypeID,)=sys.argv[1:]
phenotypeID=int(phenotypeID)-1
if(phenotypeID<0): exit("phenotype must be 1-5")

phenotypes=loadPhenotypes(PHENOTYPES)
effects=loadEffects(EFFECTS)
IDs=None
for chrom in CHROMS:
    vcf=VCF+"/"+chrom+"/"+chrom+".vcf.gz"
    print("processing",vcf,flush=True)
    IN=gzip.open(vcf,"rt")
    for line in IN:
        if(rex.find("^#CHROM",line)):
            fields=line.rstrip().split()
            IDs=fields[9:]
            continue
        if(rex.find("^\s*#",line)): continue
        fields=line.rstrip().split()
        (chrom,pos,variant,ref,alt)=fields[:5]
        genotypes=fields[9:]
        processVariant(variant,ref,alt,IDs,genotypes,phenotypes,effects,
                       phenotypeID)
    IN.close()


