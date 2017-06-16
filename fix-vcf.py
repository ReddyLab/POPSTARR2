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

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=2):
    exit(ProgramName.get()+" <in.vcf.gz>\n")
(infile,)=sys.argv[1:]

variants={}
for line in gzip.open(infile):
    line=line.decode("utf-8").rstrip()
    if(len(line)>0 and line[0]=="#"):
        print(line)
        continue
    fields=line.split("\t")
    if(len(fields)>=9):
        if(len(fields[0])<3 or fields[0][:3]!="chr"):
            fields[0]="chr"+fields[0]
        if(fields[2]=="."): fields[2]=fields[0]+"@"+fields[1]
        line=""
        for field in fields[:len(fields)-1]: line+=field+"\t"
        line+=fields[len(fields)-1]
        print(line)
    else:
        print(line)
        continue

