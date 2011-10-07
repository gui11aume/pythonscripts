# -*- coding: utf-8 -*-
"""
Set of tools put here, usually for readability of the
main code.
"""

def get_closest(dict_a, dict_b, dist=None):
   """Return a dictionary with keys as in dict_a and values
   as keys of dict_b such that the value is the element of
   dict_b closest to the key.

   The argument dist allows to compute the distance between
   element of dict_a and dict_b."""

   if dist is None:
      def dist(x, y):
         """Default infinite norm."""
         try:
            return max([abs(a-b) for (a,b) in zip(x,y)])
         except TypeError:
            return float('inf')

   # Step 1: merge the dictionaries, but for this we need
   # to ascertain that keys do not overlap.
   if set(dict_a).intersection(dict_b):
      raise KeyError('keys of dict_a and dict_b overlap')
   merge = dict(dict_a.items() + dict_b.items())
   closest = {}

   # Step 2: do a forward pass on 'merge' sorted by value.
   # Set as closest the rightmost element left of a key.
   left = None
   for key in sorted(merge, key = merge.get):
      if key in dict_a:
         closest[key] = left
      else:
         left = key

   # Step 3: do a backward pass on 'merge'. Use the dist
   # function to get closest from either left or right.
   right = left
   for key in sorted(merge, key = merge.get, reverse = True):
      if key in dict_a:
         if dist(key, right) < dist(key, closest[key]):
            closest[key] = right
      else:
         right = key

   return closest
