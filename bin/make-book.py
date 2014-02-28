import sys

# Header required to make this a Jekyll file.
HEADER = '''---
layout: book
root: .
---'''

def main(filenames):
    '''Splice a bunch of files together to make a book.'''
    print HEADER
    for f in filenames:
        with open(f, 'r') as reader:
            lines = reader.readlines()
            lines = remove_jekyll_header(lines)
            lines = remove_toc(lines)
            sys.stdout.writelines(lines)

def remove_jekyll_header(lines):
    '''Remove a dash-delimited Jekyll header (if present).'''
    if not lines[0].startswith('---'):
        return lines

    count = 0
    for (i, ln) in enumerate(lines):
        if ln.startswith('---'):
            count += 1
        if count == 2:
            break
    return lines[i+1:]

def remove_toc(lines):
    '''Remove a div with class "toc" (if present).'''
    start = end = None
    for (i, ln) in enumerate(lines):
        if ('<div class="toc"' in ln):
            start = i
        if (start is not None) and ('</div>' in ln):
            end = i+1
            break
    if (end is not None):
        lines = lines[0:start] + lines[end:]
    return lines

if __name__ == '__main__':
    main(sys.argv[1:])
