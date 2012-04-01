# The Shell

[Back To The Menu](http://github.com/thehackerwithin/UofCSCBC2012/)
- [Forward to Python Variables](http://github.com/thehackerwithin/UofCSCBC2012/tree/master/2a-PythonVariables/)

* * * * *

**Presented By : Milad Fatenejad**

**Most of this material came from a presentation by: ??????**

# What is the shell how do I access the shell?

The *shell* is a program that presents a command line interface
which allows you to control your computer using commands entered
with a keyboard instead of controlling graphical user interfaces
(GUIs) with a mouse/keyboard combination.

A *terminal* is a program you run that gives you access to the
shell. There are many different terminal programs that vary across
operating systems. To open the default terminal in Ubuntu, just
click... or just use CTRL+ALT+t
	 
There are many reasons to learn about the shell. In my opinion, the
most important reasons are that: 

1.  It is very common to encounter the shell and
    command-line-interfaces in scientific computing, so you will
    probably have to learn it eventually 

2.  The shell is a really powerful way of interacting with your
    computer. GUIs and the shell are complementary - by knowing both
    you will greatly expand the range of tasks you can accomplish with
    your computer. You will also be able to perform many tasks more
    efficiently.

The shell is just a program and there are many different shell
programs that have been developed. The most common shell (and the one
we will use) is called the Bourne-Again SHell (bash). Even if bash is
not the default shell, it usually installed on most systems and can be
started by typing "bash" in the terminal. Many commands, especially a
lot of the basic ones, work across the various shells but many things
are different. I recommend sticking with bash and learning it well.

# The Example: Manipulating Experimental Data Files

We will spend most of our time learning about the basics of the shell
by manipulating some experimental data from a hearing tests. 

**Cochlear Implants**

A cochlear implant is a small electronic device that is surgically
implanted in the inner ear to give deaf people a sense of
hearing. More than a quarter of a million people have them, but there
is still no widely-accepted benchmark to measure their effectiveness.
In order to establish a baseline for such a benchmark, our supervisor
got teenagers with CIs to listen to audio files on their computer and
report:

1.  the quietest sound they could hear
2.  the lowest and highest tones they could hear
3.  the narrowest range of frequencies they could discriminate

To participate, subjects attended our laboratory and one of our lab
techs played an audio sample, and recorded their data - when they
first heard the sound, or first heard a difference in the sound.  Each
set of test results were written out to a text file, one set per file.
Each participant has a unique subject ID, and a made-up subject name.
Each experiment has a unique experiment ID. The experiment has
collected 351 files so far.

The data is a bit of a mess! There are inconsistent file names, there
are extraneous "NOTES" files that we'd like to get rid of, and the
data is spread across many directories. We are going to use shell
commands to get this data into shape. By the end we would like to:

1.  Put all of the data into one directory called "alldata"

2.  Have all of the data files in there, and ensure that every file
    has a ".txt" extension

3.  Get rid of the extraneous "NOTES" files

If we can get through this example in the available time, we will move
onto more advanced shell topics...

# Let's get started

One very basic command is `echo`. This command is just prints text to
the terminal. Try entering the command:

    echo Hello, World

Then press enter. You should see the text "Hello, World" printed back
to you. The echo command is useful for 

## Moving around the file system

Let's learn how to move around the file system using command line
programs. This is really easy to do using a GUI (just click on
things). Once you learn the basic commands, you'll see that it is
really easy to do in the shell too. 

First we have to know where we are. The program `pwd` (print working
directory) tells you where you are sitting in the directory tree. The
command `ls` will list the files in files in the current
directory. Directories are often called "folders" because of how they
are represented in GUIs. Directories are just listings of files. They
can contain other files or directories.

[[[START AT HOME - Slide 28]]]

**File Types**

[[[FILE TYPES - Slide 29]]]

Directories can contain other files or directories. In Ubuntu, and
many operating systems, 

**Changing Directories**

[[[CHANGING DIRECTORIES - Slide 30]]]

## Command Arguments

Most programs take additional arguments that control their exact
behavior. For example, `-F` is an argument to `ls`.  The `ls` program,
like many programs, take a lot of arguments. But how do we know what
the options are to particular commands?

Most commonly used shell programs have a manual. You can access the
manual using the `man` program. Try entering:

    man ls

This will open the manual page for `ls`. Use the arrow keys to go up
and down. When you are done reading, just hit `q` to exit.



# Extra Commands

## The backtick

## Some more common commands

**which**

**alias**

## .bashrc

## ssh and scp

## Regular Expressions

# Milad's Notes:

Don't we have to clone the repo?

Introduce less early - go over searching. 
