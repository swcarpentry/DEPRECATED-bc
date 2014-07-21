import sys
import numpy as np


def main():
    script = sys.argv[0]
    action = sys.argv[1]
    filenames = sys.argv[2:]
    assert action in ['--min', '--mean', '--max', '--sum'], \
           'Action is not one of --min, --mean, or --max: ' + action
    if len(filenames) == 0:
        process(sys.stdin, action)
    else:
        for f in filenames:
            process(f, action)


def process(filename, action):
    data = np.loadtxt(filename, delimiter=',')

    data_shape = data.shape
    if len(data_shape) == 1:
        data = data.reshape((1, len(data)))

    if action == '--min':
        values = data.min(axis=1)
    elif action == '--mean':
        values = data.mean(axis=1)
    elif action == '--max':
        values = data.max(axis=1)
    elif action == '--sum':
        values = data.sum(axis=1)

    for m in values:
        print m

main()
