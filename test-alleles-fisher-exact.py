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

def getCounts(filename,variants,MIN_COUNT):
    counts={}
    with open(filename,"rt") as IN:
        for line in IN:
            fields=line.rstrip().split()
            if(len(fields)!=7): continue
            (id,chr,pos,ref,alt,refCount,altCount)=fields
            refCount=int(refCount); altCount=int(altCount)
            if(refCount+altCount<MIN_COUNT): continue
            counts[id]=[refCount,altCount]
    return counts

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=8):
    exit(ProgramName.get()+" <in.vcf.gz> <dna.counts> <rna.counts> <alpha> <all|sig> <min-count> <pseudocount>\n")
(vcf,dnaFile,rnaFile,ALPHA,allOrSig,MIN_COUNT,PSEUDOCOUNT)=sys.argv[1:]
ALPHA=float(ALPHA)
MIN_COUNT=int(MIN_COUNT)
PSEUDOCOUNT=int(PSEUDOCOUNT)

# Process the VCF file to get the ref and alt alleles
variants={}
for line in gzip.open(vcf):
    if(len(line)>0 and line[0]=="#"): continue
    fields=line.decode("utf-8").rstrip().split()
    if(len(fields)<9): continue
    (chr,pos,id,ref,alt,x,Pass,flags,GT)=fields[:9]
    altFields=alt.split(",")
    alt=altFields[0] # Keeping only the first alternate allele
    #if(id=="."): continue
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
    table=[[dnaRef,dnaAlt],[rnaRef,rnaAlt]]
    (oddsRatio,P)=stats.fisher_exact(table) # 2-sided test!
    pvalues.append(P)
    tests.append([variant,P,dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt])
(reject,q)=multipletests(pvalues,ALPHA,"fdr_bh")[:2]
print("chrom\tpos\tvariant\tP\tPadj\teffect\tDNAref\tDNAalt\tRNAref\tRNAalt\tref\talt")
for i in range(len(q)):
    if(allOrSig=="all" or (allOrSig=="sig" and q[i]<=ALPHA)):
        (variant,P,dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt)=tests[i]
        (ref,alt,chr,pos)=variants[variant]
        if(dnaRef==0 or dnaAlt==0 or rnaRef==0 or rnaAlt==0): continue
        #effectSize=(rnaAlt/dnaAlt)/(rnaRef/dnaRef) if dnaRef>0 and \
        #    dnaAlt>0 and rnaRef>0 else 0
        dnaRef+=PSEUDOCOUNT; dnaAlt+=PSEUDOCOUNT;
        rnaRef+=PSEUDOCOUNT; rnaAlt+=PSEUDOCOUNT
        effectSize=(rnaAlt/dnaAlt)/(rnaRef/dnaRef)
        print(chr,pos,variant,P,q[i],effectSize,dnaRef,dnaAlt,rnaRef,
              rnaAlt,ref,alt,sep="\t")

