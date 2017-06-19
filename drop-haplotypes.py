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
from FastaReader import FastaReader
from FastaWriter import FastaWriter
from Rex import Rex
rex=Rex()

#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=3):
    exit(ProgramName.get()+" <in.fasta> <out.fasta>\n")
(infile,outfile)=sys.argv[1:]

OUT=open(outfile,"wt")
writer=FastaWriter()
reader=FastaReader(infile)
while(True):
    (defline,seq)=reader.nextSequence()
    if(not defline): break
    if(not rex.find(">chr",defline)): continue
    writer.addToFasta(defline,seq,OUT)
OUT.close()
