#!/usr/bin/env python
"""
A C function that segfaults, and a wrapper that can handle it
"""

import ctypes

def cause_segmentation_fault():
    """Crashes the Python interpreter by segfaulting."""
    i = ctypes.c_char('a')
    j = ctypes.pointer(i)
    c = 0
    while True:
        j[c] = 'a'
        c += 1
    return j

if __name__ == "__main__":
    import sys
    
    # This handles the segfault and gives us a traceback
    if "--handle" in sys.argv:
        import faulthandler
        faulthandler.enable()
    
    # This runs the function that will segfault
    cause_segmentation_fault()