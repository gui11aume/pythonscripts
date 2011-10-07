#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import with_statement

import re
import sys

from vtrack import vheader, vskip

canonID = {}

# Read-in the lookup table in a dict.
with open(sys.argv[1]) as lookup:
   for line in vskip(lookup):
      canonID.update([line.rstrip().split('\t')])

def to_canonID(FBmatch):
   return canonID.get(FBmatch.group(), 'no_canonID')

# Write the vheader.
sys.stdout.write(vheader(*sys.argv))
# Read-in arg file and update FBgn line by line.
with open(sys.argv[2]) as argfile:
   for line in argfile:
      sys.stdout.write(
            re.sub('FBgn[0-9]{7}', to_canonID, line)
         )
