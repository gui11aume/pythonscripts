#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Write a matrix of gene counts given two .json files
containing the GO specifications and the genes in
different categories (typically targets of a protein).

The output file can be uploaded in R by using
'read.delim(file=..., comment.char="#")'.
"""

import sys
try:
   import json
except:
   import simplejson as json

from vtrack import vheader, vskip

__author__ = 'Guillaume Filion'
__email__ = 'guillaume.filion@gmail.com'
__date__ = '2011-10-10'
__version__ = '0.1'


def count(json_tar, json_specs, out=sys.stdout):
   """Print a count matrix of the target genes in a given
   GO category."""

   specs = json.load(vskip(open(json_specs)))
   tar = json.load(vskip(open(json_tar)))

   # First thing: remove the genes tagged "NA" if present.
   NA = tar.pop('NA', None)
   # There must be a "total" key of the target dictionary
   # that will be used for the counts.
   total = tar.pop('total')

   # Turn to sets to remove potential duplicate entries.
   for key in specs: specs[key] = set(specs[key])
   for key in tar: tar[key] = set(tar[key])

   # Get all GO-annotated genes.
   annotated = set([gene for ls in specs.values() for gene in ls])

   result = [] # List of count lists.

   # Proteins in line, GO terms in column. Fill by line.
   for prot in sorted(tar):
      # Count the intersection of GO and targets.
      result.append(
         [prot] + \
         [str(len(specs[GO].intersection(tar[prot]))) for GO in specs] + \
         [str(len(annotated.intersection(tar[prot])))]
      )
   # The last line contains the total genes with a GO term.
   result.append(
      ['total'] + \
      [str(len(specs[GO].intersection(total))) for GO in specs] + \
      [str(len(annotated.intersection(total)))]
   )

   # Print the table header.
   out.write('\t'.join(specs) + '\ttotal\n')
   # Print the table
   for line in result:
      out.write('\t'.join(line) + '\n')
   


if __name__ == '__main__':
   # Print the vheader.
   sys.stdout.write(vheader(*sys.argv))
   count(sys.argv[1], sys.argv[2])
