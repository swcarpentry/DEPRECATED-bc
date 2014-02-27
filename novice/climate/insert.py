import sys

assert len(sys.argv) == 3, 'Usage: insert.py country what'
country = sys.argv[1]
what = sys.argv[2]

template = "insert into readings values('{0}','{1}',{2},{3});"
for line in sys.stdin:
    year, value = line.strip().split(',')
    print template.format(country, what, year, value)
