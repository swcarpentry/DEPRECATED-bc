---
layout: lesson
root: ../..
---
# Using the IPython infrastructure for parallel computing

## Motivation

As a scientist, you are continually processing data.  This processing can
benefit greatly from parallelization, but due to various barriers and/or
habits, these techniques are often underused. For example, suppose you
are performing an experiment to determine differential expression of genes in
two different tissue types. You do all of the lab work, creating cDNA libraries,
transforming your cells, and picking your clones.  You send your clones off for
Sanger sequencing and get back from the sequencing center thousands
of chromatogram files. You now have two choices for analysis: serial or parallel.

## Objectives

After going through this lesson, you will be able to:

1. Start an IPython cluster, both manually and through the IPython notebook.
1. Connect to that cluster using the client.
1. Define which of your code you want to run in parallel.
1. Execute your code in parallel.

This lesson assumes familiarity with the Python language, navigating and
executing commands using the command line, and running Python code from
the IPython notebook.

## Introduction

The purpose of this lesson is to show you how to process your data (Sanger sequences
or otherwise) in parallel, whenever possible, and to recognize situations where
parallelization can greatly speed up and streamline your workflow. Fortunately
(or not) Python provides many packages which can meet this need, some better
suited to certain types of problems than others, with varying degrees of
user-friendliness. Some of these packages that you might have seen are:

* Symmetric multiprocessing
  * Threading
  * Multiprocessing
  * pp
* Compute Clusters
  * Celery
  * Disco
  * mpi4py

Symmetric multiprocessing (SMP) and cluster computing are two models
of parallelization running either multiple cores on a single server or
across multiple servers, respectively. Rather than having to choose
between packages, IPython lets you use either with little or no
changes to your existing code. Indeed, running your code
across multiple machines is effectively identical to running through a
local cluster from a user's perspective, outside of setting up a
remote client connected to an IPython cluster with access to engines
on remote nodes. 

## Infrastructure

IPython provides a robust and functional
infrastructure for parallel computing though multiple components
working together behind the scenes, waiting to receive and run your jobs.
These connections are managed by a [client](../../gloss.html#client) which provides connections
to a [view](../../gloss.html#view), which are the entry points to the cluster.
For the full details on the IPython parallel infrastructure,
please refer to the [official documentation](http://ipython.org/ipython-doc/stable/parallel/).

## Creating a Profile

The first thing we need to do is create a profile to store the settings
for our cluster.  This can be done by running the following command
from the command line:

~~~
$ ipython profile create --parallel --profile=cluster
~~~
{:class="in"}

This will create a directory called
`ipython_cluster` in your `IPYTHONDIR` and include configuration files for
both engine and controller in that directory. These files contain a
huge array of customization possibilities for running the cluster,
and these will be covered in subsequent lessons. For now, we will
leave these files alone since we will be working on a single machine.

## Starting the Cluster

Now that your cluster (called `cluster`) has been created, let's start it up
using its profile.  This can be done in two ways:

### Using ipcluster

From a command line, execute the `ipcluster` command to start the
cluster. Let's use a 4-node cluster:

~~~
$ ipcluster start -n 4 --profile=cluster
~~~
{:class="in"}

Note: if you run `ipcluster` without `--profile`, IPython will use
the default IPython profile.  This is fine for simple cases,
but a good habit is to always use a parallel profile,
even if you're only using default settings.

### Using the IPython Notebook

From your running IPython notebook, go to the `Clusters` tab, enter
the number of engines (i.e., nodes) in your cluster (in this case,
enter 4), and click `Start`.

## Synchronous vs. Asynchronous Tasks

The next choice is whether you process your code in
parallel or do you do it serially (synchronously). IPython lets you
do it both ways, using similar code.

### Getting the Data

To make this more concrete, imagine that you have your chromatograms in 12
different directories representing 6 different treatment conditions, each
with a replicate. The directories are labeled tX_1 and tX_2 where X is the
treatment number and 1 or 2 indicates the replicate
(i.e., the actual directory names are t1\_1, t1\_2, and so on). You've
already executed [phred](http://www.phrap.org/phredphrapconsed.html) and extracted
your sequences.  Start the IPython Notebook and execute the
following code in a cell:

~~~
!pip install wget
import wget
wget.download("http://dl.dropboxusercontent.com/u/861789/swc/lessons/ipython_parallel/Archive.zip")
~~~

Note the `!` in front of the `pip` command. This tells the IPython notebook to execute
the process in the shell (as if you were at the command line). This process should be
OS-agnostic.

Unzip the files by entering the following into a new cell. You should see the output
from the command creating the directories.

~~~
import zipfile
zipfile.ZipFile("Archive.zip").extractall()
~~~

Now that the directories are unzipped, let's store them into an array for
accessing them later, in a variable named `dirs`.

~~~
dirs = !ls | grep '^t'
~~~

We can also examine the files in each of the directories, just to have a
look at what we're dealing with:

~~~
!ls {dirs[0]}
~~~

The `{}` notation allows us to pass Python variables into the shell.

### Setting Up the Code

Now that we have access to all of the files, let's set up a function that we can
call to see how this all works. It's important to know what the distribution of
quality scores is for your sequences across positions, so let's create a function to do that:

~~~
def get_quality_distribution(seq_dir):
    import os
    positions = {}
    qual_files = [x for x in os.listdir(seq_dir) if x.endswith("qual")]
    for f in qual_files:
        f = os.path.join(seq_dir, f)
        pos = 0
        for line in open(f):
            line = line.strip()
            if line.startswith(">"):
                pos = 0
            else:
                line = [int(x) for x in line.split()]
                for elem in line:
                    if not pos in positions:
                        positions[pos] = []
                    positions[pos].append(elem)
                    pos += 1
    return positions, os.getpid()
~~~

## Coding Example

To get an idea of the time it will take to get the quality score distributions,
let's run it as one might normally, without any parallelization:

~~~
%%timeit
for d in dirs:
    get_quality_distribution(d)
~~~

Make a note of this time so that we can compare against it later.

### Create the IPython Client

Now let's add some IPython parallel magic.
In a new cell, enter the following code, and execute the cell:

~~~
from IPython.parallel import Client
rc = Client(profile="cluster")
dview = rc[:]
lview = rc.load_balanced_view()
~~~

This code sets up a client called `rc` (for remote client) using our `cluster` profile
and creates two `DirectView` instances to all of the nodes in the cluster called
`dview` and `lview`, which will enable some basic load balancing and single-host
submission of tasks.

Go back up to your `get_quality_distribution` function and add the following
`lview` annotation to the function.

~~~
import os
@lview.remote()
def get_quality_distribution(seq_dir):
    ...
~~~

This enables each of the remote engines to know about the function and allow
for its remote execution.

### Submit Your Jobs

Let's get the quality distributions for each of our 6 treatments and replicates:

~~~
dists = []
for d in dirs:
    dists.append(get_quality_distribution(d))
~~~

That was pretty fast, but is the speedup because the code was parallelized?
To find out, let's look at the types of the elements in `dists`. In a
new cell, run this:

~~~
[x for x in dists]
~~~

As the output shows, the elements of `dists` are not actually the results, but rather
`AsyncResult` objects that represent the eventual results. To get the *actual* results,
we need to use their `get` method.  Because its output is so big, we'll
just display the length of the results:

~~~
[len(x.get()) for x in dists]
~~~

Is this the result you expected? Take a look at what is returned from `get_quality_distribution`.

### Verify Remote Engines

We did something tricky, and very useful it turns out, in the return from our function: we
included the OS-level process identifier (PID). To verify that your job ran remotely, you can view all of them
from the `AsyncResult` object:

~~~
[x.get()[1] for x in dists]
~~~

If you wanted to look at the raw data, you could do `x.get()[0]` and pull the entire
dictionary, but it's just too much data to look at by hand. Let's plot it, instead.

~~~
%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

def plot_quality_scores(pos_data, title):
    plot_mean = []
    keys = sorted(pos_data.keys())
    for key in keys:
        plot_mean.append(np.mean(pos_data[key]))
    plt.plot(plot_mean)
    plt.title("%s [%d, %d]" % (title, min(plot_mean), max(plot_mean)))
    plt.xlabel("Position")
    plt.ylabel("Quality")
    plt.show()

for i, x in enumerate(dists):
    plot_quality_scores(x.get()[0], dirs[i])
~~~

## Wrap Up

To process all of my files, one at a time, would take about 3
seconds. No one would argue that this is a long time. However, when we
look at parallelizing, the jobs each individually only take 1/3 of a
second.  On 4 CPUs, we can expect all jobs to complete in about (0.333
* (12/4)) or about 1 second.  As things scale up and more processors
are added and tasks become more time consuming these benefits can grow
from seconds to days.
