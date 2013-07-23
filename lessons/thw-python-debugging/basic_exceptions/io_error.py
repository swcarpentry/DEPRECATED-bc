#!/usr/bin/env python
"""
IOError - trying to open a file that doesn't exist (usually).
"""

# This will error unless the code is in the same directory as a file named
# "a_file.txt"
with open("a_file.txt") as in_file:
    print in_file.read()
