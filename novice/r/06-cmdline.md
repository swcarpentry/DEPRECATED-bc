---
layout: lesson
root: ../..
---



## Command-Line Programs

The R Console and other interactive tools like RStudio are great for prototyping code and exploring data, but sooner or later we will want to use our program in a pipeline or run it in a shell script to process thousands of data files.
In order to do that, we need to make our programs work like other Unix command-line tools.
For example, we may want a program that reads a data set and prints the average inflammation per patient:

~~~
$ Rscript readings.R --mean inflammation-01.csv
5.45
5.425
6.1
...
6.4
7.05
5.9
~~~

but we might also want to look at the minimum of the first four lines

~~~
$ head -4 inflammation-01.csv | Rscript readings.R --min
~~~

or the maximum inflammations in several files one after another:

~~~
$ Rscript readings.py --max inflammation-*.csv
~~~

Our overall requirements are:

1. If no filename is given on the command line, read data from [standard input](../../gloss.html#standard-input).
2. If one or more filenames are given, read data from them and report statistics for each file separately.
3. Use the `--min`, `--mean`, or `--max` flag to determine what statistic to print.

To make this work, we need to know how to handle command-line arguments in a program, and how to get at standard input.
We'll tackle these questions in turn below.

<div class="objectives" markdown="1">
#### Objectives

*   Use the values of command-line arguments in a program.
*   Handle flags and files separately in a command-line program.
*   Read data from standard input in a program so that it can be used in a pipeline.
</div>

### Command-Line Arguments

Using the text editor of your choice, save the following line of code in a text file called `session-info.R `:


<div class='out'><pre class='out'><code>sessionInfo()
</code></pre></div>

The function, `sessionInfo`, outputs the version of R you are running as well as the type of computer you are using (as well as the versions of the packages that have been loaded).
This is very useful information to include when asking others for help with your R code.

Now we can run the code in the file we created from the Unix Shell using `Rscript`:


<pre class='in'><code>Rscript session-info.R</code></pre>




<div class='out'><pre class='out'><code>R version 3.1.1 (2014-07-10)
Platform: x86_64-pc-linux-gnu (64-bit)

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
 [9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  base     
</code></pre></div>

> **Tip:** If that did not work, remember that you must be in the correct directory.
You can determine which directory you are currently in using `pwd` and change to a different directory using `cd`.
For a review, see this [lesson](../shell/01-filedir.html) or the [Unix Shell Reference](../ref/01-shell.html).

Now let's create another script that does something more interesting. Write the following lines in a file named `print-args.R`:


<div class='out'><pre class='out'><code>args <- commandArgs()
args
</code></pre></div>

The function `commandArgs` extracts all the command line arguments and returns them as a vector.
Let's see what happens when we run this program in the Unix Shell:


<pre class='in'><code>Rscript print-args.R</code></pre>




<div class='out'><pre class='out'><code>[1] "/usr/lib/R/bin/exec/R" "--slave"               "--no-restore"         
[4] "--file=print-args.R"   "--args"               
</code></pre></div>

From this output, we learn that `Rscript` is just a convenience command for running R scripts.
The first argument in the vector is the path to the `R` executable.
The following are all command-line arguments that affect the behavior of R.
From the R help file:

*  `--slave`: Make R run as quietly as possible
*  `--no-restore`:  Don't restore anything that was created during the R session
*  `--file`: Run this file
*  `--args`: Pass these argments to the file being run

Thus running a file with Rscript is an easier way to run the following:


<pre class='in'><code>R --slave --no-restore --file=print-args.R --args</code></pre>




<div class='out'><pre class='out'><code>[1] "/usr/lib/R/bin/exec/R" "--slave"               "--no-restore"         
[4] "--file=print-args.R"   "--args"               
</code></pre></div>

If we run it with a few arguments, however:


<pre class='in'><code>Rscript print-args.R first second third</code></pre>




<div class='out'><pre class='out'><code>[1] "/usr/lib/R/bin/exec/R" "--slave"               "--no-restore"         
[4] "--file=print-args.R"   "--args"                "first"                
[7] "second"                "third"                
</code></pre></div>

then `commandArgs` adds each of those arguments to the vector it returns.
Since the first elements of the vector are always the same, we can tell `commandArgs` to only return the arguments that come after `--args`.
Let's update `print-args.R` and save it as `print-args-trailing.R`:


<div class='out'><pre class='out'><code>args <- commandArgs(trailingOnly = TRUE)
args
</code></pre></div>

And then run `print-args-trailing` from the Unix Shell:


<pre class='in'><code>Rscript print-args-trailing.R first second third</code></pre>




<div class='out'><pre class='out'><code>[1] "first"  "second" "third" 
</code></pre></div>

Now `commandArgs` returns only the arguments that we listed after `print-args-trailing.R`.
