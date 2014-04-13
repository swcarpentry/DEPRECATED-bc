The Shell
=========

Introduction
------------

We're here to teach you about computers and computers really do four things:

* run programs
* store data 
* communicate with each other
* and interact with us

They can do the last of these in many different ways.
The most common way of doing this is using a graphical interface.
Another way to do this is using a text-based command line interface

To interact with our computers using the command line we use something called a
command shell, which is a program whose job is to run other programs. We're
going to be using the Bash shell, which is the default shell on Mac and
Linux. Bash stands for Bourne Again SHell because it was written by Stephen
Bourne.

The shell is very useful to scientists for three main reasons:

1. Automating data analysis and data management
2. Combining existing tools in powerful ways with only a few keystrokes
3. Interacting with remote machines

I'll give you two quick examples of how we regularly use the shell in my group:

1. Who here has some aspect of their research that involves a bunch of different
   files, or a bunch of different tabs in an spreadsheet that all hold the same
   kind of information and need to be processed in similar ways? We often have a
   number of different data files that represent data from a number of different
   studies, or sites, or times. We typically need to run the same analysis on
   all of these different files and that analysis often involves combining code
   written in both Python and R, and sometimes using programs written by other
   people in other languages (e.g., like MARK for mark-recapture analyses). The
   things we're going to learn this morning make it easy to do this in a fully
   automated way, so we don't have to spend days or even weeks repeatedly
   running analyses by hand.
2. Who has done some data analysis or simulation modeling recently that took a
   long time to run and your computer was horribly slow while it was running? We
   also do simulation modeling, which requires running the same basic code with
   different sets of parameters and then combining the results. Since these
   simulations take a long time to run we typically run them on a university
   cluster so that our desktops and laptops aren't grindingly slow for days
   while the analysis runs. The things we'll learn this morning make automating
   this process and running it on someone elses computer much easier.

Some of the things we'll learn this morning can be accomplished in R or Python,
but in many cases this requires understanding how to do them in the shell.

**Open a terminal;
Windows open GitBash;
Mac type 'terminal' into the Finder;
STICKIES**

Files and Directories
---------------------

### Viewing the file system

Now that we have a shell open let's see where we are.

    pwd

This stands for "print working directory" and tells us where we are in the file
system. Your terminal probably started in your home directory like mine
did. Specifically this output tells us that I'm in a directory named `ethan`
that is a subdirectory of a directory named `home`. Windows users will probably
see a letter `c` at the beginning.

When we run `pwd` it works by:

1. finding a program called `pwd`,
2. running that program,
3. and displaying that program's output,
4. and then displaying a new prompt to let us know that it's ready for more input

To take a look at what's in the current directory use

    ls

which stands for "listing". This will print out a list of files and directories
just like using the graphical file browser on your computer and we can interact
with files and directories/folder in the usual ways using the command line.

In some cases the directories will be colored but if not we can see which names
are directories and which are files by using the -F argument

    ls -F
	
which adds a trailing slash to the directories.  Dash is how we pass options to
commands, so this says "run ls in such a way that it shows me the folders".


### Changing directories

To change directories

    cd Desktop
	
Since `Desktop` is in our current directory we can just type it's name.

If we now run

    `pwd`

and

    `ls`

We see that we are in the Desktop folder and that it is the same as our
computer's desktop. In other words, you're Desktop is just a special folder in
your home directory.

We can also specify the full path to allow us to jump anywhere.

    cd /home/ethan/
	
Windows users need to start full paths with the drive letter.

    cd /c/Users/username

If you want to move up one directory we can use

    cd ..

### Example

To help us learn the shell we're going to use an example that captures some of
the common situations we run into when analyzing data.

Nelle is a marine biologist. She recently returned from a survey collecting data
on gelatinous marine life. She has 250 samples and has assayed them all to
determine the relative abundance of 300 different proteins. Each assay results
in a single file with one line for each protein.

To analyze these data she needs to:
1. Calculate statistics for each of the proteins separately using a program her
   supervisor wrote called `goostat`.
2. Compare the statistics for each protein with corresponding statistics for
   each other protein using a program one of the other graduate students wrote
   called `goodiff`.

Running the `goostat` program 250 times would be time consuming, but
possible. Running the nearly 60K pairwise comparisons one at a time would take
at least several hundred hours.

**Let's practice what we've learned so far by looking at Nelle's data. Download
  and unzip the file on your Desktop. Use the shell to look at what is in the
  north-pacific-gyre directory and any subdirectories that it contains.**

**What did you find?**

Great. And if you don't want to have to `cd` into a directory just to see what
is there you can always use `ls` with the directory name.

    ls 2012-07-03


### Naming

Now Nelle's done a really nice job here with something. What do you think it is?

That's right, she's used names that contain a lot of information. And this is

**Rule #1. Always use meaningful, consistently structured, names.**.

Her naming conventions do two things. First they help chunk information for
us. By looking at the names of the directories we know what they
contain. And in fact each of the file names contains her labs unique code for
the location, time, depth, and other characteristics of the sample.

Second, using consistently structured file names makes it easier automatically
process the data, because it will be easy to identify the data files we want.


#### Autocomplete / Up & Down Arrows

These informative names are great, but does it seem to anyone else like they
require a lot of typing.

    cd /h<tab>/e<tab>/D<tab>
	
If nothing happens then either there is nothing to complete or there is too
much. A second ``<tab>`` will show you the available options.

You can also use the up and down arrows to cycle through commands and edit only
the part of a command that you want to change.

**Rule #2. If the computer can do it, let the computer do it.**


**Challenge?**


### Making directories

We can add new directories using

    mkdir thesis

Running

    ls

shows us that thesis exists,

    ls thesis

but that there's nothing in it.
	

### Making files

To make files in the shell we use simple text editors.
I'm going to use `nano` today since it's one of the most basic and since it should
work for everyone.

To make a new file we simply type ``nano`` and then the name of the file.
Let's make a rough outline of our thesis:

    nano outline.txt

Type in a few lines of text:

```
1. Collect data
2. ?
3. Profit
```

then use Control-O to save our data (the O stands for out, the ^ in most Unix
documentation means Control).

Control-X quits the editor and return to the shell. `nano` doesn't leave any
output on the screen after it exits, but

    ls

now shows that we have created a file called `outline.txt`.
	
You can also use your favorite text editor for this, just make sure to save your
file in the right directory.

*Demonstrate the same thing using gedit creating a new data2.txt file*

**RED LIGHT/GREEN LIGHT**

### Copying, Moving, Renaming, and Deleting Files

We can make copies of files using ``cp`` for copy

    cp outline.txt copy_of_outline.txt

where first parameter tells `cp` what we're copying, and the second is where to
put the copy and what to call it.
	
	ls

shows us that we've made a copy of the file in the thesis directory.

And we can move files using ``mv``

    mkdir temp
	ls
	mv copy_of_outline.txt temp
	ls
	cd temp
	ls

We can also use mv to rename files using ``mv`` by "moving" them to a new
file name. Maybe we're feeling enthusiastic about the state of our outline so we
decide to rename it.

    mv outline.txt rough_draft.txt
	ls

If we want to get rid of files we use `rm` for remove

    rm rough_draft.txt


**Create a backup copy of your thesis, you do keep backups of your data right?;
Make a new outline.txt file;
Make a directory called backup;
Make a copy of your outline.txt file;
Move it into the backup folder;
Change the backup file's name to include todays date**



**BREAK: When we come back we'll learn about how to combine different programs**


### Shell Programs and Redirects

Now that we know how to work with the shell we can use it to do powerful things
by combining lots of small built in programs.

Nelle's assays have now finished running, so she can get back to working with
the data itself.

The first thing she wants to do is run some basic sanity checks on the
data. Remember that each assay is supposed to include measures of 300 different
proteins, one per line, so Nelle wants to look to see if any files have too many
or too few lines.

She can find out the number of lines in a file using the command `wc` for "word count".

    wc NENE01729A.txt

This actually shows us the number of characters, words, and lines in the file,
so we can add a `-l` for "lines" to tell `wc` to just show the line counts.

    wc -l NENE01729A.txt

That's great for one file, but we need the counts for all the files. `wc` will
take multiple file names as input

    wc -l NENE01729A.txt NENE01729B.txt

but that's a lot of typing even with the tab key.


### Wildcards

We can use wildcards to work with multiple files at once. So

    wc -l *.txt

counts the lines for all of the files ending in .txt. The shell expands the
wildcard before executing `wc`, so this is identical to typing out all of the
filenames.

With a small number of files we could probably just look through the counts to
see if there were any problems, but with hundreds or thousands of files that
won't work very well, so let's store this output so that we can work with it
further.

    wc -l *.txt > lengths

The `>` sign redirects the output of the `wc` command into a file called
lengths, instead of printing it to the screen.

We can look at the contents of `lengths` using

    cat lengths

`cat` stands for "concatenate": it prints the contents of files one after
another. There's only one file in this case, so cat just shows us what it
contains.

Now let's use the `sort` command to sort its contents.

    sort lengths

And we can save that to a file called `sorted-lengths`

    sort lengths > sorted-lengths

Finally we can look at just the last few lines of our sorted file to see if
there are any files with too many lines.

    tail sorted-lengths

**Perform the same basic santity check, but only for the files that end in A.txt
and use a command called head, which does the opposite of tail, to show the
files with the fewest rows**

**What did you see? Any problems**

When Nelle goes back and checks it, she sees that she did that assay at 8:00 on a
Monday morning—someone was probably in using the machine on the weekend, and she
forgot to reset it.

### Pipes

We can also skip the intermediate files by using pipes.
Pipes pass the data that would have gone into the intermediate program
straight to the next program.

So, the commands we just wrote:

    sort lengths > sorted-lengths
    tail sorted-lengths

can be written like this:

    sort lengths | tail

The output from the first command is passed through the pipe to the second
command which takes it was input. This simple idea is why Unix has been so
successful. Instead of creating enormous programs that try to do many different
things, Unix programmers focus on creating lots of simple tools that each do one
job well, and that work well with each other. This programming model is called
pipes and filters. We've already seen pipes; a filter is a program like wc or
sort that transforms a stream of input into a stream of output.


### Other people's programs

Piping together unix tools can be very powerful for automating data analysis.  I
often use this for quick sanity checks on much more complicated Python code. But
the real power of the piping in my every day life is piping my tools and other
scientists tools together. This makes it easy to use tools written in languages
you're not familiar with and to combine tools from different languages.

**Rule #3. Don't reinvent the wheel.**

So, recall that Nelle needs to calculate some statistics for each of her
datasets. She could figure out how to do this herself, but one of her lab mates
already wrote code to do this for another project. Unfortunately the code is
written in Python and Nelle is only familiar with R, but that's OK, because she
can make it work using the shell.

The file is called `goostats.py`. To run this code we need to tell the shell to
run it using python, which we do by giving it the name of the program that will
run it, then the name of our program, and then the input.

    python goostats.py NENE01729A.txt
	
This can then be integrated into our pipeline. So if we want to sort based
on the total number of individuals:
	
    python goostats.py NENE01729A.txt | head -1
	

**Using pipes run goostats on only the data from the first 50 rows of
  NENE01729A.txt and store the results to a file with a name indicating the
  analysis that the file contains. The make a directory called analysis and move
  the file into that directory**

**BREAK**


### Loops

This is great for a single datafile but Nelle has a lot of datafiles and needs
to analyze them all at once. `goostats` is only designed to handle one file at a
time, so to do this we need to use a for loop. This loop will do some operation
once for each thing in a list. In this case the list of things is a list of
filenames, so we'll start by figuring out how to just print that list of file
names.

    for datafile in *.txt
	do
	    echo $datafile
	done

When the shell sees the keyword `for`, it knows it is supposed to repeat a command
(or group of commands) once for each thing in a list. In this case, the list is
the two filenames. Each time through the loop, the name of the thing currently
being operated on is assigned to the variable called filename. Inside the loop,
we get the variable's value by putting `$` in front of it: `$filename` is
NENE01729A.txt the first time through the loop, NENE01729B.txt the second, and so
on. Finally, the command that's actually being run tells the shell to print, so
this loop prints out the name of each data file in turn.

There's nothing special about the word datafile...

Our for loop can include as many commands as we want

    for datafile in *.txt
	do
	    echo $datafile
		head $datafile
	done


**COPY TO ETHERPAD
Write a for loop that, for each datafile, prints the name of the file and then
runs goostats on it**

Really want to store the output

    for datafile in *.txt
	do
        python goostats.py $datafile > stats-$datafile
	done

And so now we have all of the output of our statistics runs in a well named and
structured format that is ready for further analysis.


### Bash scripting

Once we have the commands that we want working we typically want to store them
since if we needed to do it once, we'll probably need to do it again.  To do
this easily we can store them in a shell script. We do this by adding the
commands to a text file and then running that text file from the command line.

**Who remembers all of the commands we used exactly to get to this point?**

History

    for datafile in *.txt
	do
        python goostats.py $datafile > stats-$datafile
	done

    bash do-stats.sh

Now, as Nelle finishes up more assays she can simply add the datafiles to her
directory and rerun this function and automatically update her analysis. In
addition, if Nelle or her advisor ever change their minds about the details of
the analysis, they can make a simple change and then rerun everything.

**What's missing?** (docs)

**Rule #4 Document the purpose of your code** This will make it much easier for
  other people (including your future self) to re-use the code

**Modify your shell script so that the output only includes the first line of
output from goostats, the all important 'goo-coefficient'**

Let's add one more thing to our script before we're done. Often we don't want to
store our raw data and the output from our analyses in the same directory.

So let's design our bash script so that it can take an argument that tells it where
we want to store the output.

To do this we use a special variable `$1`. Inside a shell script, $1 means "the
first filename (or other parameter) on the command line".
	
    mkdir $1
    for datafile in *.txt
	do
        python goostats.py $datafile > $1/stats-$datafile
	done

 We can now run our script like this:

    bash do-stats.sh goostats-results

We've built this script up piece by piece, and in fact that's exactly how most
people work and in fact I find that it is much more efficient to do something
small, make sure that it's working, and then do the next part, then to try to
do everything at once and then spend hours tracking down all of the bugs.
