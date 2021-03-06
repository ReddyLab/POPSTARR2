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
from Pipe import Pipe
import gzip
from Rex import Rex
rex=Rex()
from scipy import stats
from statsmodels.stats.multitest import multipletests

MAX_VARIANTS=-1
BETA_BINOMIAL="/home/bmajoros/src/util/beta-binomial"
USE_FISHER=False

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

def betaBinomial(dnaRef,dnaAlt,rnaRef,rnaAlt):
    OPTIONS=" -t -s 2 "
    cmd=BETA_BINOMIAL+""+OPTIONS+" "+str(rnaAlt)+" "+str(rnaRef+rnaAlt)+" "+\
        str(dnaAlt+1)+" "+str(dnaAlt+dnaRef+2)
    #cmd=BETA_BINOMIAL+""+OPTIONS+" "+str(rnaAlt)+" "+str(rnaRef+rnaAlt)+" "+\
    #    str(dnaRef+1)+" "+str(dnaAlt+dnaRef+2)
    #print(cmd)
    output=Pipe.run(cmd)
    #print("OUTPUT:",output)
    P=float(output)
    #print(P)
    return P

def getCounts(filename,variants,MIN_COUNT):
    counts={}
    numLoaded=0
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=7): continue
            (id,chr,pos,ref,alt,refCount,altCount)=fields
            refCount=int(refCount); altCount=int(altCount)
            if(refCount+altCount<MIN_COUNT): continue
            counts[id]=[refCount,altCount]
            numLoaded+=1
            if(MAX_VARIANTS>0 and numLoaded>=MAX_VARIANTS): break
    return counts

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=9):
    exit(ProgramName.get()+" <in.vcf.gz> <dna.counts> <rna.counts> <alpha> <all|sig> <min-count> <pseudocount> <want-power:0|1>\n")
(vcf,dnaFile,rnaFile,ALPHA,allOrSig,MIN_COUNT,PSEUDOCOUNT,WANT_POWER)=\
    sys.argv[1:]
ALPHA=float(ALPHA)
MIN_COUNT=int(MIN_COUNT)
PSEUDOCOUNT=int(PSEUDOCOUNT)
WANT_POWER=int(WANT_POWER)

# Process the VCF file to get the ref and alt alleles
variants={}
for line in gzip.open(vcf):
    if(len(line)>0 and line[0]=="#"): continue
    fields=line.decode("utf-8").rstrip().split()
    if(len(fields)<9): continue
    (chr,pos,id,ref,alt,x,Pass,flags,GT)=fields[:9]
    altFields=alt.split(",")
    alt=altFields[0] # Keeping only the first alternate allele
    if(id=="."): id=chr+"@"+pos
    variants[id]=[ref,alt,chr,pos]

# Process DNA and RNA files
rnaCounts=getCounts(rnaFile,variants,MIN_COUNT)
dnaCounts=getCounts(dnaFile,variants,MIN_COUNT)

# Test for differences
pvalues=[]
tests=[]
for variant in dnaCounts.keys():
    (ref,alt,chr,pos)=variants[variant]
    dnaRec=dnaCounts[variant]
    rnaRec=rnaCounts.get(variant,None)
    if(rnaRec is None): continue
    (dnaRef,dnaAlt)=dnaRec
    (rnaRef,rnaAlt)=rnaRec
    P=None
    if(USE_FISHER):
        table=[[dnaRef,dnaAlt],[rnaRef,rnaAlt]]
        (oddsRatio,P)=stats.fisher_exact(table) # 2-sided test!
    else: P=betaBinomial(dnaRef,dnaAlt,rnaRef,rnaAlt)
    power=None
    if(WANT_POWER):
        power=estimatePower(dnaRef,dnaAlt,rna,foldChange,iterations)
    pvalues.append(P)
    tests.append([variant,P,dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt,power])
(reject,q)=multipletests(pvalues,ALPHA,"fdr_bh")[:2]
print("chrom\tpos\tvariant\tpower\tP\tPadj\teffect\tDNAref\tDNAalt\tRNAref\tRNAalt\tref\talt")
for i in range(len(q)):
    if(allOrSig=="all" or (allOrSig=="sig" and q[i]<=ALPHA)):
        (variant,P,dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt,power)=tests[i]
        (ref,alt,chr,pos)=variants[variant]
        if(dnaRef==0 or dnaAlt==0 or rnaRef==0 or rnaAlt==0): continue
        dnaRef+=PSEUDOCOUNT; dnaAlt+=PSEUDOCOUNT;
        rnaRef+=PSEUDOCOUNT; rnaAlt+=PSEUDOCOUNT
        effectSize=(rnaAlt/dnaAlt)/(rnaRef/dnaRef)
        powerString="NA" if power is None else str(power)
        print(chr,pos,variant,powerString,P,q[i],effectSize,dnaRef,dnaAlt,
              rnaRef,rnaAlt,ref,alt,sep="\t")

