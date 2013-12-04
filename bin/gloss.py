#!/usr/bin/env python

'''
Check glossary entries.
Usage: gloss.py glossary_file html_file...
Typically, bin/gloss.py ./gloss.md _site/*/novice/*.html
'''

import sys
import re

def main(gloss_filename, html_filenames):
    '''Main driver.'''
    known = get_gloss_entries(gloss_filename)
    for f in html_filenames:
        report_missing(f, known)
    report_unused(known)

def get_gloss_entries(filename):
    '''Get entries from glossary, reporting any that are out of order or
    duplicated.  Result is a dictionary of anchors to counts
    (initially all 0).  Checks along the way that internal definitions
    resolve.'''

    # Entry pattern: 0 = key, 1 = abbrev, 2 = term
    p_entry = re.compile(r'\*\*([^\*]+)\*\*(\s+\(.+\))?:\s+<a\s+name="([^"]+)"></a>')

    # Use pattern: 0 = key
    p_use = re.compile(r'\([^\)]+\)\[\#([^\]]+)\]')

    result = {}
    internal = set()
    undone = 0
    last_seen = ''
    out_of_order = []

    with open(filename, 'r') as reader:
        for line in reader:
            m = p_entry.search(line)
            if m:
                text = m.group(1)
                if text == 'FIXME':
                    undone += 1
                else:
                    if text.lower() < last_seen:
                        out_of_order.append(text)
                    last_seen = text.lower()
                key = m.group(3)
                if key in result:
                    print 'Duplicate key {0} in {1}'.format(key, filename)
                result[key] = 0
            for ref in p_use.findall(line):
                internal.add(ref)

    if undone:
        print '{0} UNDONE'.format(undone)

    if out_of_order:
        print 'OUT OF ORDER:'
        for term in out_of_order:
            print '   ', term

    missing_internal = internal - set(result.keys())
    if missing_internal:
        print 'MISSING (INTERNAL):'
        for term in sorted(missing_internal):
            print '   ', term

    return result

def report_missing(f, known):
    '''Read HTML files to find glossary definitions not in the glossary
    file.  Counts the number of times each term used so that unused
    glossary entries can be reported.'''

    # Use pattern: 0 == upward ref to glossary file, 1 == key
    p_use = re.compile(r'<a href="(\.\./)*gloss.html#([^"]+)">')

    with open(f, 'r') as reader:
        unknown = set()
        for line in reader:
            matches = p_use.findall(line)
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
