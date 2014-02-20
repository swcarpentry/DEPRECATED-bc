#!/usr/bin/env python

import sys
while True:
    line = sys.stdin.readline()
    if not line:
        break
    if ('Entity' in line):
        sys.stdin.readline()
        sys.stdin.readline()
    else:
        sys.stdout.write(line)
