#!/usr/bin/env python
"""
Setting a trace in some sample code with pdb.
"""

import pdb

a = 2
while True:
    print "This number's gonna get HUGEE"
    pdb.set_trace()
    a *= a
