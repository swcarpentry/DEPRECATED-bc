#!/usr/bin/env python
"""
TypeError - trying to do an operation with incompatible data types.
"""

# You can add strings or ints together, no problem
added_strings = "1" + "1"
added_ints = 1 + 1

# But what is the sum of an int and a string? (spoiler: it's a TypeError)
added_mixed_types = added_strings + added_ints

# Here's something to mess with you. Does this work?
multiplied_mixed_types = added_strings * added_ints