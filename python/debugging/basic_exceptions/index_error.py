#!/usr/bin/env python
"""
IndexError - accessing a list by an index that doesn't exist. This usually
comes up due to some messed up logic around the source of the error.
"""

# Try accessing an element of this list
a_list = ["you", "should", "avoid", "accessing", 
          "lists", "by", "index", "anyway!"]
non_existent_element = a_list[100]

# Pythonic bonus points: only use a list if you're going to be iterating over
# every element. The following won't cause an IndexError.
for item in a_list:
    print item

# If you need to access random elements of a list, consider rewriting your code
# to use a dictionary.