#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the Gnu General Public
# License (GPL) version 3, as described at www.opensource.org.
# Copyright (C)2016 William H. Majoros (martiandna@gmail.com).
#=========================================================================
from __future__ import (absolute_import, division, print_function,
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)

import scipy as sp
from scipy import special
import numpy as np

#=========================================================================
# Attributes:
#   
# Instance Methods:
#   BetaBinomial()
# Class Methods:
#   
#=========================================================================
class BetaBinomial:
    """BetaBinomial"""
    def __init__(self):
        pass

    def P(this,alpha,beta,n,k):
        part_1 = sp.special.comb(n,k)
        part_2 = sp.special.betaln(k+alpha,n-k+beta)
        part_3 = sp.special.betaln(alpha,beta)
        result = (np.log(part_1) + part_2)- part_3
        return np.exp(result)

