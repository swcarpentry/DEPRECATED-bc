import fileinput

def get_raw_data():
    """Load data from stdin"""
    data = []
    for line in fileinput.input():
        data.append(float(line.strip()))
    return data

def main():
    data = get_raw_data()
    print "goo-coefficient:", sum(data) / len(data)
    print "wilsons-q:", sum([value ** 2 for value in data])
    print "omega:", len(data)
    

if __name__ == '__main__':
    main()
