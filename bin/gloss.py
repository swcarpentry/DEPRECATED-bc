#!/usr/bin/env python

'''
Check glossary entries.
Usage: gloss.py glossary_file html_file...
Typically, bin/gloss.py ./gloss.md _site/*/novice/*.html
'''

import sys
import re

# glossary entries: first group is the text, second is the anchor
GLOSS_ENTRY_PAT = re.compile(r'\*\*([^\*]+)\*\*:\s+<a\s+name="([^"]+)"></a>')

# glossary definition: first group is the key
GLOSS_USE_PAT = re.compile(r'<a href="(\.\./)*gloss.html#([^"]+)">')

def main(gloss_filename, html_filenames):
    '''Main driver.'''
    known = get_gloss_entries(gloss_filename)
    for f in html_filenames:
        report_missing(f, known)
    report_unused(known)

def get_gloss_entries(filename):
    '''Get entries from glossary, reporting any that are out of order or
    duplicated.  Result is a dictionary of anchors to counts
    (initially all 0).'''
    result = {}
    last_seen = ''
    out_of_order = []
    with open(filename, 'r') as reader:
        for line in reader:
            m = GLOSS_ENTRY_PAT.search(line)
            if m:
                text = m.group(1)
                if text < last_seen:
                    out_of_order.append(text)
                last_seen = text
                key = m.group(2)
                if key in result:
                    print 'Duplicate key {0} in {1}'.format(key, filename)
                result[key] = 0
    if out_of_order:
        print 'OUT OF ORDER:'
        for term in out_of_order:
            print '   ', term
    return result

def report_missing(f, known):
    '''Read HTML files to find glossary definitions not in the glossary
    file.  Counts the number of times each term used so that unused
    glossary entries can be reported.'''
    with open(f, 'r') as reader:
        unknown = set()
        for line in reader:
            matches = GLOSS_USE_PAT.findall(line)
            for (prefix, m) in matches:
                if m in known:
                    known[m] += 1
                else:
                    unknown.add(m)
        if unknown:
            print 'UNDEFINED FROM {0}'.format(f)
            for term in sorted(unknown):
                print '   ', term

def report_unused(known):
    '''Report unused glossary entries.'''
    temp = [k for k in known if not known[k]]
    if not temp:
        return
    print 'UNUSED'
    for term in sorted(temp):
        print '   ', term

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2:])
