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


## Saving time with shortcuts, wild cards, and tab completion

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

**Wild cards**

Navigate to the `~/UofCSCBC2012/Shell-1/data/THOMAS` directory. This
directory contains our hearing test data for THOMAS. If we type `ls`,
we will see that there are a bunch of files which are just four digit
numbers. By default, `ls` lists all of the files in a given
directory. The `*` character is a shortcut for "everything". Thus, if
you enter `ls *`, you will see all of the contents of a given
directory. Now try this command:

    ls *1

This lists every file that ends with a `1`. This command:

    ls /usr/bin/*.sh

Lists every file in `/usr/bin` that ends in the characters `.sh`. And
this command:

    ls *4*1

lists every file in the current directory which contains the number
`4`, and ends with the number `1`. There are four such files: `0241`,
`0341`, `0431`, and `0481`. 

So how does this actually work? Well...when the shell (bash) sees a
word that contains the `*` character, it automatically looks for files
that match the given pattern. In this case, it identified four such
files. Then, it replaced the `*4*1` with the list of files, separated
by spaces. In other the two commands:

    ls *4*1
    ls 0241 0341 0431 0481

are exactly identical. The `ls` command cannot tell the difference
between these two things.

* * * *
**Short Exercise**

Do each of the following using a single `ls` command without
navigating to a different directory.

1.  List all of the files in `/bin` that contain the letter `a`
2.  List all of the files in `/bin` that contain the letter `a` or the letter `b`
3.  List all of the files in `/bin` that contain the letter `a` AND the letter `b`
* * * *

**Tab Completion**

Navigate to the home directory. Typing out directory names can waste a
lot of time. When you start typing out the name of a directory, then
hit the tab key, the shell will try to fill in the rest of the
directory name. For example, enter:

    ls U<tab>

The shell will fill in the rest of the directory name for
`UofCSCBC2012`. Now enter:

    ls D<tab><tab>

When you hit the first tab, nothing happens. The reason is that there
are multiple directories in the home directory which start with
D. Thus, the shell does not know which one to fill in. When you hit
tab again, the shell will list the possible choices. 

Tab completion can also fill in the names of programs. For example,
enter `e<tab><tab>`. You will see the name of every program that
starts with an `e`. One of those is `echo`. If you enter `ec<tab>` you
will see that tab completion works.

## Which program? ##

Commands like `ls`, `rm`, `echo`, and `cd` are just ordinary programs
on the computer. A program is just a file that you can *execute*. The
program `which` tells you the location of a particular program. For
example:

    which ls

Will return "/bin/ls". Thus, we can see that `ls` is a program that
sits inside of the `/bin` directory. Now enter:

    which mplayer

You will see that `mplayer` is a program that sits inside of the
`/usr/bin` directory.

So ... when we enter a program name, like `ls`, and hit enter, how
does the shell know where to look for that program? How does it know
to run `/bin/ls` when we enter `ls`. The answer is that when we enter
a program name and hit enter, there are a few standard places that the
shell automatically looks. If it can't find the program in any of
those places, it will print an error saying "command not found". Enter
the command:

    echo $PATH

This will print out the value of the `PATH` environment variable. More
on environment variables later. Notice that a list of directories,
separated by colon characters, is listed. These are the places the
shell looks for programs to run. If your program is not in this list,
then an error is printed. The shell ONLY checks in the places listed
in the `PATH` environment variable. 

Navigate to the `1-Shell` directory and list the contents. You will
notice that there is a program (executable file) called `hello` in
this directory. Now, try to run the program by entering:

    hello

You should get an error saying that hello cannot be found. That is
because the directory `/home/thw/UofCSCBC2012/1-Shell` is not in the
`PATH`. You can run the `hello` program by entering:

    ./hello

Remember that `.` is a shortcut for the current working
directory. This tells the shell to run the `hello` program which is
located right here. So, you can run any program by entering the path
to that program. You can run `hello` equally well by specifying:

    /home/thw/UofCSCBC2012/1-Shell/hello

Or by entering:

    ../1-Shell/hello

When there are no `/` characters, the shell assumes you want to look
in one of the default places for the program.


## Examining Files

We now know how to switch directories, run programs, and look at the
contents of directories, but how do we look at the contents of files?

The easiest way to examine a file is to just print out all of the
contents using the program `cat`. Enter the following command:

    cat ex_data.txt

This prints out the contents of the `ex_data.txt` file. If you enter:

    cat ex_data.txt ex_data.txt

It will print out the contents of `ex_data.txt` twice. `cat` just
takes a list of file names and writes them out one after another (this
is where the name comes from, `cat` is short for concatenate). 

* * * *
**Short Exercises**

1.  Print out the contents of the `/usr/share/dict/american-english`
    file. What does this file contain?

2.  Without changing directories, (you should still be in `1-Shell`),
    use one short command to print the contents of all of the files in
    the /home/milad/UofCSCBC2012/1-Shell/data/THOMAS directory.
* * * *

`cat` is a terrific program, but when the file is really big, it can
be annoying to use. The program, `less`, is useful for this
case. Enter the following command:

    less /usr/share/dict/american-english

`less` opens the file, and lets you navigate through it. The commands
are identical to the `man` program. Use "space" to go forward and hit
the "b" key to go backwards. The "g" key goes to the beginning of the
file and "G" goes to the end. Finally, hit "q" to quit.

`less` also gives you a way of searching through files. Just hit the
"/" key to begin a search. Enter the name of the word you would like
to search for and hit enter. It will jump to the next location where
that word is found. Try searching the `american-english` file for the
word "copper". If you hit "/" then "enter", `less` will just repeat
the previous search. `less` searches from the current location and
works its way forward. If you are at the end of the file and search
for the word "copper", `less` will not find it. You need to go to the
beginning of the file and search.

Remember, the `man` program uses the same commands, so you can search
documentation using "/" as well!

* * * *
**Short Exercise**

Use the commands we've learned so far to figure out what the `-fs`
argument for the program `mplayer` does. `mplayer` video playing program.
 * * * * 


## Redirection

We now know a lot of the basic shell commands. Let's turn to the
experimental data from the hearing tests that we began with. This data
is located in the `~/UofCSCBC2012/1-Shell/data` directory. Each
subdirectory corresponds to a particular participant in the
study. Navigate to the `Bert` subdirectory in `data`.  There are a
bunch of text files which contain experimental data results. Lets
print them all:

    cat au*

Now enter the following command:

    cat au* > ../all_data

This tells the shell to take the output from the `cat au*` command and
dump it into a new file called `../all_data`. To verify that this
worked, examine the `all_data` file. If `all_data` had already
existed, we would overwritten it. So the `>`

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
