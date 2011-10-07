#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Write JSON gene target list for binding data.
"""

import sys
try:
   import json
except ImportError:
   import simplejson as json

# Private modules.
from vtrack import vheader
from misctools import get_closest

__author__ = 'Guillaume Filion'
__email__ = 'guillaume.filion@gmail.com'
__date__ = '2011-10-07'
__version__ = '0.1'

# Beyond that distance (in bp), a TSS has no closest
# binding element (array probe).
MAXDIST = 1000

def dist(TSS, mapinfo):
   """Distance function between a TSS and the mapping
   information of a binding element. Passed to
   'get_closest()'."""
   if TSS[0] != mapinfo[0]:
      # Different seqnames: infinite distance.
      return float('inf')
   else:
      # Same seqname: distance to the closest end of mapinfo
      return min(
            abs(TSS[1] - mapinfo[1]),
            abs(TSS[1] - mapinfo[2])
         )


def JSONtargets(mappingfile, bindingfile):
   """TODO: write a docstring."""
   # Read in gene mapping. Skip comment lines and remove stray
   # 'chr' sometimes present in chromosome names.
   mapping = [
         l.rstrip().replace('chr','').split('\t') \
         for l in open(mappingfile, 'r') \
         if l[0] != '#'
      ]

   # Remove the header if present (recognized by 'start' and
   # 'end' in third and fourth columns.
   if mapping[0][2:4] == ['start','end']: mapping.pop(0)

   # Collect TSS, if gene is on +, TSS is on start, else on end.
   TSS = {}
   for row in mapping:
      thisTSS = {
        '+': lambda x: (x[1], int(x[2])), # 2nd and 3rd column.
        '-': lambda x: (x[1], int(x[3]))  # 2nd and 4th column.
      }.get(row[4])(row)
      # Arrange geneIDs by TSS in a dictionary.
      # Example: TSS['FBgn0031208'] = ('2L', 7529)
      TSS[row[0]] = thisTSS


   # Read in binding data. Skip comment lines and remove
   # 'chr' on chromosome names.
   binding = [
         l.rstrip().replace('chr','').split('\t') \
         for l in open(bindingfile, 'r') \
         if l[0] != '#'
      ]
   # Get feature names and remove (pop) the header.
   # Example: features = ['D005', 'D007', ...]
   features = binding.pop(0)[4:]
   targets = {'NA': []}
   for feature in features:
      targets[feature] = []


   # Collect mapping information (seqname, start, end) and
   # binding info (0/1).
   mapinfo = {}
   bindinfo = {}
   for row in binding:
      # Example: mapinfo['r5GATC2L00037'] = ('2L', 5301, 6026)
      mapinfo[row[0]] = (row[1], int(row[2]), int(row[3]))
      # Example: bindinfo['r5GATC2L00037'] = [0,0,1,...]
      bindinfo[row[0]] = row[4:]


   # Get the closest feature to TSS.
   close_elt = get_closest(TSS, mapinfo, dist = dist)


   for geneID in close_elt:
      if dist(TSS[geneID], mapinfo[close_elt[geneID]]) > MAXDIST:
         # The gene is too far. Push it to NA.
         targets.get('NA').append(geneID)
      else:
         # The gene gets the status of the binding element closest
         # to its TSS.
         for feature in [
               feat for (feat, yes) in \
               # Example: [('D005', 0), ('D007', 0), ...]
               zip(features, bindinfo[close_elt[geneID]]) \
               if yes == '1'
            ]:
            targets.get(feature).append(geneID)


   # Print the version tracking header and the JSON data.
   sys.stdout.write(vheader(*sys.argv))
   json.dump(targets, sys.stdout, indent=4)


if __name__ == '__main__':
   """Call on a pair of files.
   Input #1: gene mapping table [geneID, seqname, start, end, strand]
   Input #2: binding table [ID, seqname, start, end, feature1, ...]"""
   JSONtargets(sys.argv[1], sys.argv[2])
