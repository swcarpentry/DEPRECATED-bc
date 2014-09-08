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

With this in hand, let's build a version of `readings.py` that always prints the per-patient (per-row) mean of a single data file.
The first step is to write a function that outlines our implementation, and a placeholder for the function that does the actual work.
By convention this function is usually called `main`, though we can call it whatever we want.
Write the following code in a file called `readings-01.R`:


<div class='out'><pre class='out'><code>main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  dat <- read.csv(file = filename, header = FALSE)
  mean_per_patient <- apply(dat, 1, mean)
  return(mean_per_patient)
}
</code></pre></div>


This function gets the name of the file to process from the first element returned by `commandArgs`.
Here's a simple test to run from the Unix Shell:


<pre class='in'><code>Rscript readings-01.R inflammation-01.csv</code></pre>

There is no output because we have defined a function, but haven't actually called it.
Let's add a call to `main` and save it as `readings-02.R`:


<div class='out'><pre class='out'><code>main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  filename <- args[1]
  dat <- read.csv(file = filename, header = FALSE)
  mean_per_patient <- apply(dat, 1, mean)
  return(mean_per_patient)
}

main()
</code></pre></div>


<pre class='in'><code>Rscript readings-02.R inflammation-01.csv</code></pre>




<div class='out'><pre class='out'><code> [1] 5.450 5.425 6.100 5.900 5.550 6.225 5.975 6.650 6.625 6.525 6.775 5.800
[13] 6.225 5.750 5.225 6.300 6.550 5.700 5.850 6.550 5.775 5.825 6.175 6.100
[25] 5.800 6.425 6.050 6.025 6.175 6.550 6.175 6.350 6.725 6.125 7.075 5.725
[37] 5.925 6.150 6.075 5.750 5.975 5.725 6.300 5.900 6.750 5.925 7.225 6.150
[49] 5.950 6.275 5.700 6.100 6.825 5.975 6.725 5.700 6.250 6.400 7.050 5.900
</code></pre></div>

#### Challenges

  + Write a command-line program that does addition and subtraction:


<pre class='in'><code>Rscript arith.R 1 + 2</code></pre>




<div class='out'><pre class='out'><code>[1] 3
</code></pre></div>


<pre class='in'><code>Rscript arith.R 3 - 4</code></pre>




<div class='out'><pre class='out'><code>[1] -1
</code></pre></div>



  + What goes wrong if you try to add multiplication using `*` to the program?
  


