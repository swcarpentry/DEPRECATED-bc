#!/usr/bin/env python

'''Make a temporary page displaying all the images used in a set of lessons.'''

import sys
import os
import re

HEADER = '''<html>
<body>'''

ENTRY = '<div><img src="{0}" alt="{1}" /><p>{1} ({0})</p></div>'

FOOTER = '''</body>
</html>'''

MARKDOWN_P = re.compile(r'<img\s+src="([^"]+)"\s+alt="([^"]+)"\s*/>')
IPYNB_P = re.compile(r'<img\s+src=\\"([^"]+)\\" alt=\\"([^"]+)\\"\s*/>')

def main(filenames):
    print HEADER
    for f in filenames:
        with open(f, 'r') as reader:
            for line in reader:
                display(f, line)
    print FOOTER

def display(filename, line):
    for p in (MARKDOWN_P, IPYNB_P):
        m = p.search(line)
        if not m: continue
        relative_path = m.group(1)
        alt_text = m.group(2)
        modified_path = adjust_path(filename, relative_path)
        print ENTRY.format(modified_path, alt_text)

def adjust_path(base_filename, relative_path):
    fixed = os.path.join(os.path.dirname(base_filename), relative_path)
    return fixed.replace('/files/', '/')

if __name__ == '__main__':
    main(sys.argv[1:])
