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
import gzip
from Rex import Rex
rex=Rex()
from scipy import stats
from statsmodels.stats.multitest import multipletests

def getCounts(filename,variants):
    counts={}
    IN=open(filename,"rt")
    for line in IN:
        fields=line.rstrip().split()
        if(len(fields)<2): continue
        id=fields[0]
        if(id=="."): continue
        variant=variants.get(id,None)
        if(variant is None): continue
        hash=counts.get(id,None)
        if(hash is None): hash=counts[id]={}
        (ref,alt)=variant[:2]
        for field in fields[1:]:
            if(not rex.find("(\S+)=(\d+)",field)):
                raise Exception("can't parse "+field)
            allele=rex[1]
            count=int(rex[2])
            if(allele==ref): hash["ref"]=count
            elif(allele==alt): hash["alt"]=count
    IN.close()
    return counts

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=5):
    exit(ProgramName.get()+" <in.vcf.gz> <dna.counts> <rna.counts> <alpha>\n")
(vcf,dnaFile,rnaFile,ALPHA)=sys.argv[1:]
ALPHA=float(ALPHA)

# Process the VCF file to get the ref and alt alleles
variants={}
for line in gzip.open(vcf):
    if(len(line)>0 and line[0]=="#"): continue
    fields=line.decode("utf-8").rstrip().split()
    if(len(fields)<9): continue
    (chr,pos,id,ref,alt,x,Pass,flags,GT)=fields[:9]
    if(id=="."): continue
    variants[id]=[ref,alt,chr,pos]

# Process DNA and RNA files
dnaCounts=getCounts(dnaFile,variants)
rnaCounts=getCounts(rnaFile,variants)

# Test for differences
pvalues=[]
tests=[]
for variant in dnaCounts.keys():
    (ref,alt,chr,pos)=variants[variant]
    dnaRec=dnaCounts[variant]
    rnaRec=rnaCounts.get(variant,None)
    if(rnaRec is None): continue
    dnaRef=dnaRec.get("ref",0)
    dnaAlt=dnaRec.get("alt",0)
    rnaRef=rnaRec.get("ref",0)
    rnaAlt=rnaRec.get("alt",0)
    table=[[dnaRef,dnaAlt],[rnaRef,rnaAlt]]
    (oddsRatio,P)=stats.fisher_exact(table)
    pvalues.append(P)
    tests.append([variant,P,dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt])
(reject,q)=multipletests(pvalues,ALPHA,"fdr_bh")[:2]
print("chrom\tpos\tvariant\tP\tPadj\tDNAref\tDNAalt\tRNAref\tRNAalt\tref\talt")
for i in range(len(q)):
    #if(q[i]<=ALPHA):
    if(True):
        (variant,P,dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt)=tests[i]
        (ref,alt,chr,pos)=variants[variant]
        print(chr,pos,variant,P,q[i],dnaRef,dnaAlt,rnaRef,rnaAlt,ref,alt,sep="\t")

