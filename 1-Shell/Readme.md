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
operating systems.
	 
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
started by typing `bash` in the terminal. Many commands, especially a
lot of the basic ones, work across the various shells but many things
are different. I recommend sticking with bash and learning it well.

To open a terminal, just double click on the "Konsole" icon on the
Desktop.

# The Example: Manipulating Experimental Data Files

We will spend most of our time learning about the basics of the shell
by manipulating some experimental data from a hearing tests. To get
the data for this test, you will need internet access. Just enter the
command:

    git clone git://github.com/thehackerwithin/UofCSCBC2012.git

This will grab all of the data needed for this workshop from the
internet.

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

Whenever you start up a terminal, you will start in a special
directory called the *home* directory. Every user has their own home
directory where they have full access to do whatever they want. In
this case, the `pwd` command tells us that we are in the `/home/thw`
directory. This is the home directory for the `thw` user. That is our
user name. You can always find out your user name by entering the
command `whoami`. 

**File Types**

When you enter the `ls` command lists the contents of the current
directory. There are several items in the home directory, notice that
they are all colored blue. This tells us that all of these items are
directories as opposed to files.

Lets create an empty file using the `touch` command. Enter the
command:

    touch testfile

Then list the contents of the directory again. You should see that a
new entry, called `testfile`, exists. It is colored white meaning that
it is a file, as opposed to a directory. The `touch` command just
creates an empty file. 

Some terminals will not color the directory entries in this very
convenient way. In those terminals, use `ls -F` instead of `ls`. The
`-F` argument modifies the results so that a slash is placed at the
end of directories. If the file is *executable* meaning that it can be
run like a program, then a star fill be placed of the file name.

You can also use the command `ls -l` to see whether items in a
directory are files or directories. `ls -l` gives a lot more
information too, such as the size of the file and information about
the owner. If the entry is a directory, then the first letter will be
a "d". The fifth column shows you the size of the entries in
bytes. Notice that `testfile` has a size of zero.

Now, let's get rid of `testfile`. To remove a file, just enter the
command:

    rm testfile

The `rm` command can be used to remove files. If you enter `ls` again,
you will see that `testfile` is gone.


**Changing Directories**

Now, let's move to a different directory. The command `cd` (change
directory) is used to move around. Let's move into the `UofCSCBC2012`
directory. Enter the following command:

    cd UofCSCBC2012

Now use the `ls` command to see what is inside this directory. You
will see that there is an entry which is green. This means that this
is an executable. If you use `ls -F` you will see that this file ends
with a star.

This directory contains all of the material for this boot camp. Now
move to the directory containing the data for the shell tutorial:

    cd 1-Shell

If you enter the `cd` command by itself, you will return to the home
directory. Try this, and then navigate back to the `1-Shell`
directory.

## Arguments

Most programs take additional arguments that control their exact
behavior. For example, `-F` is an argument to `ls`.  The `ls` program,
like many programs, take a lot of arguments. But how do we know what
the options are to particular commands?

Most commonly used shell programs have a manual. You can access the
manual using the `man` program. Try entering:

    man ls

This will open the manual page for `ls`. Use the arrow keys to go up
and down. When you are done reading, just hit `q` to exit.

Programs that are run from the shell can get extremely complicated. To
see an example, open up the manual page for the `gcc` program, which
is an open source C compiler and see just how complicated programs can
get! No one can possibly learn all of these arguments, of course. So
you will probably find yourself referring back to the manual page
frequently.


## Examining Files ##

The easiest way to examine a file ...

* * * *
**Short Exercise**

Use the commands we've learned so far to figure out what the `-Wall`
argument for `gcc` does.
* * * *

# Extra Commands

## The backtick, xargs

## Some more common commands

**which**

**alias**

**touch**

**du**

## .bashrc

## ssh and scp

## Regular Expressions

# Milad's Notes:

Don't we have to clone the repo?

Introduce less early - go over searching. 

# Background, Foreground, control-Z, control-C

## Not everything is a file or a directory...
- Symbolic links
- /dev

## Permissions

## Variables