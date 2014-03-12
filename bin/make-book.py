from __future__ import print_function
import sys
import os.path

# Header required to make this a Jekyll file.
HEADER = '''---
layout: book
root: .
---'''

def main():
    print(HEADER)
    for filename in sys.argv[1:]:
        with open(filename, 'r') as reader:
            lines = reader.readlines()

        title = None
        if lines[0].startswith('---'):
            lines, skipped = skip(filename, lines, '---', '---')
            title = extract_title(filename, skipped)

        lines, _ = skip(filename, lines, '<div class="toc"', '</div>')

        lines = fix_image_paths(filename, lines)
        lines = fix_gloss(filename, lines)

        if title:
            print(format_title(filename, title))
        for line in lines:
            print(line.rstrip())

        print()

def skip(filename, lines, open, close):
    '''Skip a block of lines starting with open and ending with close.'''
    i_open = None
    i_close = None
    for (i, ln) in enumerate(lines):
        if (i_open is None) and ln.startswith(open):
            i_open = i
        elif (i_open is not None) and ln.startswith(close):
            i_close = i
            return lines[:i_open] + lines[i_close+1:], lines[i_open:i_close]
    else:
        return lines, None

def fix_image_paths(filename, lines):
    '''Modify image paths to include directory.'''
    front, _ = os.path.split(filename)
    front = front.replace('cached/', '')
    src = '<img src="'
    dst = '<img src="{0}/'.format(front)
    for (i, ln) in enumerate(lines):
        lines[i] = ln.replace(src, dst)
    return lines

def fix_gloss(filename, lines):
    '''Fix up glossary entries.'''
    is_glossary = 'gloss.md' in filename
    for (i, ln) in enumerate(lines):
        lines[i] = ln.replace('href="../../gloss.html#', 'href="#g:')
        if is_glossary:
            lines[i] = ln.replace('](#', '](#g:').replace('<a name="', '<a name="g:')
    return lines

def extract_title(filename, lines):
    '''Extract title from YAML header.'''
    for ln in lines:
        if ln.startswith('title:'):
            return ln.split(':', 1)[1].strip()
    return None

def format_title(filename, title):
    title = '## {0}'.format(title)
    f = os.path.split(filename)[-1]
    if f in ('index.md', 'intro.md'):
        return '\n'.join(['<div class="chapter" markdown="1">', title, '</div>'])
    else:
        return title

if __name__ == '__main__':
    main()
