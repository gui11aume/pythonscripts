#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Print a table of gene mapping information.
Input: FlyBase fasta file (dmel-all-genes...)
Output: [geneID, seqname, start, end, strand]
"""

import sys
import re
import datetime
import hashlib
from vtrack import vheader

# Manual version control.
from git import *

__version__ = '0.1'
__author__ = 'Guillaume Filion'
__email__ = 'guillaume.filion@gmail.com'

# Automatic version control.
filehead = """# Script: %s
# Version: %s
# Date: %s
# Invocation: %s
# MD5: %s
"""

# Patterns.
loc = 'loc=([^:]+):(complement\()?([0-9]+)\.\.([0-9]+)'
ID = 'ID=([^;]+);'
# Get the headers of fasta file (start with '>').
f = open(sys.argv[1], 'r')
headers = [line.rstrip() for line in f if line[0] == '>']

genes = {}

for h in headers:
   (seq,strand,start,end) = re.search(loc, h).groups()
   strand = '+' if strand is None else '-'
   (geneID,) = re.search(ID, h).groups()
   # Format the line for output, update dict 'genes'.
   genes[(seq, int(start))] = '\t'.join((geneID,seq,start,end,strand))

code = open(__file__, 'r').read()
sys.stdout.write(filehead % (
      __file__,
      version??,
      datetime.datetime().strfime('%Y-%m-%d'),
      ' '.join(sys.argv),
      hashlib.md5(code).hexdigest()
   )

# Print the vheader.
sys.stdout.write(vheader(__file__))
# Sort lines by key, ie seqname, start.
for key in sorted(genes):
   sys.stdout.write('%s\n' % genes[key])
