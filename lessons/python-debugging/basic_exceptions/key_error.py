#!/usr/bin/env python
"""
KeyError - accessing a dictionary by a key that doesn't exist. Like a name
error, this is usually due to a typo.
"""

# Define a dictionary and then access a nonexistent element
a_dict = {
          "a really long key name": 1, 
          "some other": "stuff not relevant",
          "to the": "error"
         }
value = a_dict["a raelly long key name"]
