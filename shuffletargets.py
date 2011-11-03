#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
   import json
except ImportError:
   import simplejson as json
import sys
import random
from vtrack import vheader, vskip

def shuffle(jsontar):
   """Modify jsontar in place, using the 'total' field."""
   sample = random.Random().sample
   for prot in jsontar:
      jsontar[prot] = sample(jsontar['total'], len(jsontar[prot]))

if __name__ == '__main__':
   jsontar = json.load(vskip(open(sys.argv[1])))
   shuffle(jsontar)
   sys.stdout.write(vheader(*sys.argv))
   json.dump(jsontar, sys.stdout, indent=4)
