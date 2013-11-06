import sys

count = 0
for line in sys.stdin:
    count += 1

print '{0} lines in standard input'.format(count)
