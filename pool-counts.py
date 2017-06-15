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
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
counts={}
for line in sys.stdin:
    fields=line.rstrip().split()
    id=fields[0]
    alleles=counts[id]
    if(alleles is None): alleles=counts[id]={}
    for field in fields[1:]:
        if(not rex.find("(\S+)=(\d+)",field)):
            raise Exception("can't parse field: "+field)
        allele=rex[1]
        count=int(rex[2])
        alleles[allele]=alleles.get(allele,0)+count
for variant in counts.keys():
    print(variant+"\t",end="")
    alleles=counts[variant]
    for allele in alleles.keys():
        count=alleles[allele]
        print("\t"+allele+"="+str(count))
    print()

