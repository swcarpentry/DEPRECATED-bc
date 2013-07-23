#!/usr/bin/env python
"""
Examples of things pyflakes will pick up
"""

# Unused imports
from math import sin

# Using undefined variables
# (usually checking undefined variables solves typo issues)
x = no_one_defined_me + 1

# Using an uninitialized variable
non_existent_dict["a_field"] = 100
