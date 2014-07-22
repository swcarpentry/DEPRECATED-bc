#!/usr/bin/env python

import os, sys, errno
import re
import argparse
from time import time
import multiprocessing

import numpy as np
import matplotlib.pyplot as plt

def plotData(outputDir, plotNum):
	outFilename = "plot_%d.pdf" % (plotNum,)
	outFilepath = os.path.join(outputDir, outFilename)
	
	# Plot some random data
	# Adapted from: http://matplotlib.org/examples/shapes_and_collections/scatter_demo.html
	N = 500
	# First we need to re-initialize the random number generator for each worker
	# See: https://groups.google.com/forum/#!topic/briansupport/9ErDidIBBFM
	np.random.seed( int( time() ) + plotNum )
	x = np.random.rand(N)
	y = np.random.rand(N)
	area = np.pi * (15 * np.random.rand(N))**2 # 0 to 15 point radiuses

	print("\tMaking plot %d" % (plotNum,) )
	plt.scatter(x, y, s=area, alpha=0.5)
	plt.savefig(outFilepath)
	# Clear figure so that the next plot this worker makes will not contain
	# data from previous plots
	plt.clf() 
	
	return (plotNum, outFilepath)


if __name__ == '__main__':
    # Handle command line options
    parser = argparse.ArgumentParser(description='Plot random data in parallel')
    parser.add_argument('-o', '--outputDir', required=True,
                        help='The directory to which plot files should be saved')
    parser.add_argument('-n', '--numPlots', required=False, type=int, default=32,
    					help='The number of plots to make')
    parser.add_argument('--numProcessors', required=False, type=int, 
    					default=multiprocessing.cpu_count(),
    					help='Number of processors to use. ' + \
    					"Default for this machine is %d" % (multiprocessing.cpu_count(),) )
    args = parser.parse_args()

    if not os.path.isdir(args.outputDir) or not os.access(args.outputDir, os.W_OK):
    	sys.exit("Unable to write to output directory %s" % (args.outputDir,) )
    
    if args.numPlots < 1:
    	sys.exit('Number of plots must be greater than 0')
    
    if args.numProcessors < 1:
    	sys.exit('Number of processors to use must be greater than 0')
    
    # Start my pool
    pool = multiprocessing.Pool( args.numProcessors )

    print("Making %d plots of random data using %d processors..." % \
    		(args.numPlots, args.numProcessors) )

    # Build task list
    tasks = []
    plotNum = 0
    while plotNum < args.numPlots:
    	plotNum += 1
    	tasks.append( (args.outputDir, plotNum, ) )
    
    # Run tasks
    results = [pool.apply_async( plotData, t ) for t in tasks]

    # Process results
    for result in results:
        (plotNum, plotFilename) = result.get()
        print("Result: plot %d written to %s" % (plotNum, plotFilename) )

    pool.close()
    pool.join()
