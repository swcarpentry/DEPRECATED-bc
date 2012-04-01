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
behavior. For example, `-F` and `-l` are arguments to `ls`.  The `ls`
program, like many programs, take a lot of arguments. But how do we
know what the options are to particular commands?

Most commonly used shell programs have a manual. You can access the
manual using the `man` program. Try entering:

    man ls

This will open the manual page for `ls`. Use the space key to go
forward and b to go backwards. When you are done reading, just hit `q`
to exit.

Programs that are run from the shell can get extremely complicated. To
see an example, open up the manual page for the `mplayer` program,
which is command line driven video player. There are about 300
arguments to the mplayer command. No one can possibly learn all of
these arguments, of course. So you will probably find yourself
referring back to the manual page frequently.

**Examining the contents of other directories**

By default, the `ls` commands lists the contents of the working
directory (i.e. the directory you are in). You can always find the
directory you are in using the `pwd` command. However, you can also
give `ls` the names of other directories to view. Navigate to the
home directory if you are not already there. Then enter the
command:

    ls UofCSCBC2012

This will list the contents of the `UofCSCBC2012` directory without
you having to navigate there. Now enter:

    ls UofCSCBC2012/1-Shell

This prints the contents of `1-Shell`. The `cd` command works in a
similar way. Try entering:

    cd UofCSCBC2012/1-Shell

and you will jump directly to `1-Shell` without having to go through
the intermediate directory.

## Full vs. Relative Paths

The `cd` command takes an argument which is the directory
name. Directories can be specified using either a *relative* path a
full *path*. The directories on the computer are arranged into a
hierarchy. The full path tells you where a directory is in that
hierarchy. Navigate to the home directory. Now, enter the `pwd`
command and you should see:

    /home/thw

which is the full name of your home directory. This tells you that you
are in a directory called `thw`, which sits inside a directory called
`home` which sits inside the very top directory in the hierarchy. The
very top of the hierarchy is a directory called `/` which is usually
referred to as the *root directory*. So, to summarize: `thw` is a
directory in `home` which is a directory in `/`.

Now enter the following command:

    cd /home/thw/UofCSCBC2012/1-Shell

This jumps to `1-Shell`. Now go back to the home directory. We saw
earlier that the command:

    cd UofCSCBC2012/1-Shell

had the same effect - it took us to the `1-Shell` directory. But,
instead of specifying the full path
(`/home/thw/UofCSCBC2012/1-Shell`), we specified a *relative path*. In
other words, we specified the path relative to our current
directory. A full path always starts with a `/`. A relative path does
not. You can usually use either a full path or a relative path
depending on what is most convenient. If we are in the home directory,
it is more convenient to just enter the relative path since it
involves less typing.

Now, list the contents of the /bin directory. Do you see anything
familiar in there?


## Saving time with shortcuts and wild cards

**Shortcuts**
There are some shortcuts which you should know about. Dealing with the
home directory is very common. So, in the shell the tilde character,
`~`, is a shortcut for your home directory. Navigate to the `1-Shell`
directory, then enter the command:

    ls ~

This prints the contents of your home directory, without you having to
type the full path. The shortcut `..` always refers to the directory
above your current directory. Thus: 

    ls ..

prints the contents of the `/home/thw/UofCSCBC2012`. You can chain
these together, so:

    ls ../../

prints the contents of `/home/thw` which is your home
directory. Finally, the special directory `.` always refers to your
current directory. So, `ls`, `ls .`, and `ls ././././.` all do the
same thing, they print the contents of the current directory. This may
seem like a useless shortcut right now, but we'll see when it is
needed in a little while.

To summarize, the commands `ls ~`, `ls ~/.`, `ls ../../`, and `ls
/home/thw` all do exactly the same thing. These shortcuts are not
necessary, they are provided for your convenience.



## Examining Files

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