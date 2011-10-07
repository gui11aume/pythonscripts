#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Make a lookup table of FBgn synonyms from a FlyBase
dmel-all-gene... fasta file.

"""

import sys
import re
from vtrack import vheader


__version__ = '0.1'
__author__ = 'Guillaume Filion'
__email__ = 'guillaume.filion@gmail.com'



# Patterns.
syn = 'FBgn[0-9]{7}'
ID = 'ID=([^;]+);'
# Get the headers of fasta file (start with '>').
f = open(sys.argv[1], 'r')
headers = [line.rstrip() for line in f if line[0] == '>']

names = {}

for h in headers:
   # Get the canonical ID.
   (geneID,) = re.search(ID, h).groups()
   # Update the geneID/canonID dict (note that we have
   # a canonID/canonID term).
   names.update(dict.fromkeys(re.findall(syn, h), geneID))

lookup_header = '\t'.join((
      'geneID',
      'canonID'
   )) + '\n'

# Print the vheader.
sys.stdout.write(vheader(__file__))
# Print a lookup header.
sys.stdout.write(lookup_header)
# Print the dictionary.
for key in sorted(names):
   sys.stdout.write('%s\t%s\n' % (key, names[key]))
