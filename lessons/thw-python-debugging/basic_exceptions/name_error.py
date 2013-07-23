#!/usr/bin/env python
"""
NameError - usually happens because of typos. Can also happen by forgetting
to initialize certain variable types.
"""

# A typo-induced NameError
a_really_complicated_variable_name = 1
something_else = a_really_complciated_name + 1

# A NameError from forgetting to initialize a list (or dictionary) data type
for i in range(100):
    my_dict[4 * i] = i - 1