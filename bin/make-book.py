import sys

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
            lines = remove(lines, '---', '---')
            lines = remove(lines, '<div class="toc"', '</div>')
            sys.stdout.writelines(lines)

def remove(lines, opener, closer):
    '''Remove lines that lie between two markers.'''
    keeping = True
    result = []
    for ln in lines:
        if keeping:
            if opener in ln:
                keeping = False
            else:
                result.append(ln)
        else:
            if closer in ln:
                keeping = True
    return result

if __name__ == '__main__':
    main(sys.argv[1:])
