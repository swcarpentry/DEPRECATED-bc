import sys
import os.path

def main():
    for filename in sys.argv[1:]:
        with open(filename, 'r') as reader:
            lines = reader.readlines()

        title = extract_title(filename, lines)
        lines = skip(filename, lines, True, '---', '---')
        lines = skip(filename, lines, False, '<div class="toc"', '</div>')

        if title:
            print make_title(filename, title)
        for line in lines:
            print line.rstrip()

        print

def extract_title(filename, lines):
    for ln in lines:
        if ln.startswith('title:'):
            return ln.split(':', 1)[1].strip()
    return None

def skip(filename, lines, required, open, close):
    i_open = None
    i_close = None
    for (i, ln) in enumerate(lines):
        if (i_open is None) and ln.startswith(open):
            i_open = i
        elif (i_open is not None) and ln.startswith(close):
            i_close = i
            return lines[:i_open] + lines[i_close+1:]
    assert not required, 'Did not find "{0}" to "{1}" in {2}'.format(open, close, filename)
    return lines

def make_title(filename, title):
    title = '## {0}'.format(title)
    f = os.path.split(filename)[-1]
    if f in ('index.md', 'intro.md'):
        return '\n'.join(['<div class="chapter">', title, '</div>'])
    else:
        return title

if __name__ == '__main__':
    main()
