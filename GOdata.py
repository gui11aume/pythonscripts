#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Provides basic class and functions for handling GO data.

Call from the command line to produce a JSON file with GO
specifications, ie an indexed list of gene ID for each GO
term. The script can be run on the (decompacted) targets of
the following links:

http://flybase.org/static_pages/downloads/...gene_association.fb.gz
http://archive.geneontology.org/latest-termdb/go_daily-termdb.obo-xml.gz
[Note that the FlyBase link name changes every release].

The script has not been tested on association files from
different species, but it should work after minor modifications.
The header of the FlyBase association files starts with '!',
which is the default 'comment_char' argument of the function
'parseGeneAssociations'. The function assumes that the gene ID
is in file column 2  and the GO term in column 5.
"""

import re
import sys
from optparse import OptionParser
from xml.sax import handler, make_parser

try:
   import json
except ImportError:
   # Python 2.5
   import simplejson as json

from vtrack import vheader, vskip


__author__ = 'Guillaume Filion'
__email__ = 'guillaume.filion@gmail.com'
__date__ = '2011-10-17'
__version__ = '0.1'




#####################################################
####                  CLASSES                    ####
#####################################################



class GOException(Exception):
   pass


class GOspecs:
   """Contains OBO-type and gene association data.
   Provides tree upward traversal.

   Assumes that namespaces are:
     cellular_component (GO:0005575)
     molecular_function (GO:0003674)
     biological_process (GO:0008150)

   Get the subslim, bu 'self.childrenOf(self.slim)'.
   To know the ancestor of GO:xxx in the slim, call
   'self.parentsIn(GO:xxx, self.slim)."""

   ###############   DATA    ################

   namespaces = ['GO:0005575', 'GO:0003674', 'GO:0008150']
   namespaceSets = {
         'GO:0005575': set(),
         'GO:0003674': set(),
         'GO:0008150': set()
      }


   
   ############## CONSTRUCTOR  ##############

   def __init__(self,
         associations, # [('gene 1', 'GO 1'), ('gene 2', 'GO 2'), ... ]
         parentDict,   # {'GO 1': ['GO 2', 'GO 3', 'GO 4], ... }
         canonid={},   # {'GO 1': 'GO 2', 'GO 3': 'GO 4, ... }
         names={},     # {'GO 1': 'name 1', ... }
         slim=[],
         assoVersion={},
         OBOversion={}
      ):
      """Initializes the following attributes:
      parentDict: child GO-term/parent GO-term dict.
      associations: (GO-term, gene) list.
      names: GO-term/name dict.
      specs: GO-term/gene list dict.
      slim: GO-term list.
      namespaceSets: namespace/GO-term list dict.
      version: string list."""

      # Store the original data.
      self.parentDict = parentDict
      self.names = names
      # Make sure the namespaces are in parentDict and
      # have no parent.
      for GOterm in self.namespaces:
         self.parentDict[GOterm] = []
      # Cast associations in-place (with canonical GO terms).
      self.associations = [list(a) for a in associations]
      for pair in self.associations:
         pair[0] = canonid.get(pair[0], pair[0])
      # Check that GO terms are in the OBO tree.
      if set([
               GO for (GO,gen) in self.associations
            ]).difference(parentDict):
         # Ooops...
         raise GOException(set([
               GO for (GO,gen) in self.associations
            ]).difference(parentDict))
     
      # Versions and dates.
      self.version = []
      self.version.extend(assoVersion)
      self.version.extend(OBOversion)

      # Build specifications.
      self.specs = {}
      for (GOterm, gene) in self.associations:
         for ancestor in self.traverseDepthFirst(GOterm):
            try:
               self.specs[ancestor].add(gene)
            except KeyError:
               self.specs[ancestor] = set([gene])
            # Update 'namespaceSets' in the process.
            if ancestor in self.namespaces:
               # namespaceSets = {'GO:0005575': {gene1, ... }, ... }
               self.namespaceSets[ancestor].add(GOterm)

      # Remove multiple entries by list-setting.
      for GOterm in self.specs:
         self.specs[GOterm] = list(set(self.specs[GOterm]))

      # Remove the name spaces from the slim by convention.
      # This allows to check that a GO term has no slim
      # ancestor. Also remove the slim entries without
      # entry in specs (e.g. no photosynthesis in Drosophila).
      self.slim = list(
            set(slim).intersection(self.specs).difference(self.namespaces)
         )


   ###########   TREE METHODS   ###########

   def traverseBreadthFirst(self, GOterm):
      """Traverse the parentDict tree breadth-first 'upward'
      from a given GOterm."""

      yield GOterm
      last = GOterm
      for node in self.traverseBreadthFirst(GOterm):
         for parent in self.parentDict[node]:
            yield parent
            last = parent
         # Stop iterations when recursion stalls.
         if last == node:
            return

   def traverseDepthFirst(self, GOterm):
      """Traverse the parentDict tree depth-first 'upward'
      from a given GOterm."""

      yield GOterm
      for parent in self.parentDict[GOterm]:
         for node in self.traverseDepthFirst(parent):
            yield node

   def childrenOf(self, GOlist):
      """Return a list of the (direct) children GO terms of
      a GO term set."""
      
      # Allow query by single string.
      GOlist = set([GOlist]) if type(GOlist) is str else set(GOlist)
      return list(set([
         GOterm for GOterm in self.specs \
         if GOlist.intersection(self.parentDict[GOterm])
      ]))

   def parentsIn(self, GOterm, termset):
      """Return the parents of a GO term in a given set of
      GO terms."""

      return list(set([
            ancestor \
            for ancestor in self.traverseDepthFirst(GOterm)
         ]).intersect(termset))



#####################################################
####                  PARSERS                    ####
#####################################################


def parseOBOXML(filename):
   """Return a couple of dictionaries and lists upon
   parsing an OBOXML file with a SAX parser."""

   # Define a SAX parser.
   class OBOXMLHandler(handler.ContentHandler):
      """SAX parser for OBO-XML file. Assumes that every "term"
      node has only one "id" subnode, and that the "is_obsolete"
      term comes last."""

      def __init__(self, OBOversion, parentDict, canonid, names, slim):
         """Pass in hashes to be updated while parsing."""
         self.id = None
         self.name = None
         self.read_to = False
         self.in_source = False
         self.data = ''
         self.OBOversion = OBOversion
         self.parentDict = parentDict
         self.canonid = canonid
         self.names = names
         self.slim = slim

      def startElement(self, name, attrs):
         self.data = ''
         if self.in_source:
            self.name = name
            return
         if name == 'source':
            self.in_source = True

      def endElement(self, name):
         if name == 'relationship':
            self.read_to = False
         elif name == 'source':
            self.in_source = False

         def no_action(): return
         def _id(): self.id = self.data.rstrip()
         def names(): self.names.setdefault(self.id, self.data.rstrip())
         def alt_id(): self.canonid[self.data.rstrip()] = self.id
         # Obsolete terms are kept, but the connection
         # to parents is severed.
         def is_obsolete(): self.parentDict[self.id] = []
         def subset():
            if self.data == 'goslim_generic':
               self.slim.append(self.id)
         def is_a():
            self.parentDict[self.id] = \
                  self.parentDict.get(self.id, []) + [self.data.rstrip()]
         def _type():
            if self.data == 'part_of': self.read_to = True
         def to():
            if self.read_to:
               is_a()

         # Action depends on the tag.
         {

            'id': _id,
            'name': names,
            'alt_id': alt_id,
            'is_obsolete': is_obsolete,
            'subset': subset,
            'is_a': is_a,
            'type': _type,
            'to': to

         }.get(name, no_action)()

         if self.in_source and self.name:
            self.OBOversion.append('%s: %s' % (name, self.data.rstrip()))

      def characters(self, data):
         # SAX stream can cut anywhere. Process when tag closes.
         self.data += data

   # End of SAX parser definition.


   OBOversion = ['-- OBO version information --']
   parentDict = {}
   canonid = {}
   names = {}
   slim = []

   parser = make_parser()
   parser.setContentHandler(OBOXMLHandler(
         OBOversion,
         parentDict,
         canonid,
         names,
         slim
   ))
   parser.parse(open(filename))

   return {
         'parentDict': parentDict,
         'canonid': canonid,
         'names': names,
         'slim': slim,
         'OBOversion': OBOversion
      }



def parseGeneAssociations(filename, comment_char='!', columns=(2,5)):
   assoVersion = ['-- associations version information --']
   pairlist = []
   for line in vskip(open(filename)):
      if line.startswith(comment_char):
         assoVersion.append(line[1:].rstrip())
      else:
         # Parse by specified columns.
         items = line.split('\t')
         gene = items[columns[0]-1]
         # Skip gene with no canonical ID (flanked by '__')
         # and lines with the 'NOT' keyword
         if gene[:2] == '__' or re.search('NOT', line):
            continue
         GOterm = items[columns[1]-1]
         pairlist.append((GOterm, gene))

   return {
         'associations': pairlist,
         'assoVersion': assoVersion
      }



#####################################################
####                    MAIN                     ####
#####################################################

def buildGO(association_fname, OBOXML_fname):
   """Wrapper to instantiate GOspecs from two files."""

   # Parse Gene associations and OBOXML.
   kwargs = parseGeneAssociations(association_fname)
   kwargs.update(parseOBOXML(OBOXML_fname))
   # Create the GOspecs instance.
   return GOspecs(**kwargs)

if __name__ == '__main__':
   """Read in two files and dump the specs to JSON."""
   
   # Parse options.
   parser = OptionParser(
         usage = '%prog [--slim] [--subslim] association_file OBOXML_file',
         description = 'Basic tool to build GO associations.',
         version = '%prog ' + __version__
      )
   parser.add_option(
         '-s',
         '--slim',
         action = 'store_true',
         default = False,
         help = 'Output GO slim specifications.'
      )
   parser.add_option(
         '--subslim',
         action = 'store_true',
         default = False,
         help = 'Output GO sub+slim specifications.'
      )
   parser.add_option(
         '-n',
         '--namespace',
         dest = 'namespace',
         default = None,
         help =
         """Restrict to a specific GO namespace. One of:
           GO:0005575 (cellular_component)
           GO:0003674 (molecular_function)
           GO:0008150 (biological_process)"""
      )

   (options, args) = parser.parse_args()

   if options.slim and options.subslim:
      sys.exit('specify either --slim or --subslim')
   if not options.namespace in \
         (None, 'GO:0005575', 'GO:0003674', 'GO:0008150'):
      sys.exit('namespace must be GO:0005575, GO:0003674 or GO:0008150')
   
   # Instantiate GOspecs.
   data = buildGO(args[0], args[1])

   # Select what to dump.
   if options.slim:
      # User specified --slim.
      specs = dict([(k, data.specs[k]) for k in data.slim])
   elif options.subslim:
      # User specified --subslim.
      subslim = set(data.childrenOf(data.slim) + data.slim)
      specs = dict([(k, data.specs[k]) for k in subslim])
   else:
      specs = data.specs
   # Restrict to a namespace if specified.
   if options.namespace:
      key_subset = set(specs).intersection(
            data.namespaceSets[options.namespace])
      specs = dict([(k, specs[k]) for k in key_subset])

   # Write the vheader.
   sys.stdout.write(vheader(*sys.argv))
   # Write the association/OBO header.
   for info in data.version:
      sys.stdout.write('# %s\n' % info)
   # Dump.
   json.dump(specs, sys.stdout, indent=4)
