#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Replace FBgn IDs in a file according to a lookup table.
If the FBgn is not in the table, it is flanked by '__'
which is something to grep for.
"""

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
   return canonID.get(
         FBmatch.group(),
         '__' + FBmatch.group()
      )

# Write the vheader.
sys.stdout.write(vheader(*sys.argv))
# Read-in arg file and update FBgn line by line.
with open(sys.argv[2]) as argfile:
   for line in argfile:
      sys.stdout.write(
            re.sub('FBgn[0-9]{7}', to_canonID, line)
         )
