def main():
    script = sys.argv[0]
    filenames = sys.argv[1:]
    for f in filenames:
        process(f)

def process(filename):
    print filename
