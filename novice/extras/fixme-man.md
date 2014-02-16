---
layout: lesson
root: ../..
title: Manual Pages
level: novice
---
You can get help for any Unix command with the `man` (short for manual)
command.  As an example here is the command to lookup information on 'cp'

~~~
$ man cp
~~~

The output will be shown directly in the terminal. The content displayed is
sometimes refered to as the "man page" (for example, you can say *there is an
option to `cp` to create links instead of copies, but I cannot remember what it
is called; you should check the man page*).

The way in which the man page is displayed depends on your system and
configuration. On a modern system, you should be able to use the Up and Down
keys to navigate (on an older system, you will only be able to see a page at at
time). When you have finished reading you can return to the command prompt with
'q'.

The output of the man command is typically very complete, but concise, as it is
designed to be used as a reference. It does not provide tutorial-like
introductions to the Unix system. For quick browsing, it is divided into
sections:

*   NAME - 
    gives the name of the command and a brief description
*   SYNOPSIS - 
    how to run the command, including optional and mandatory information.
    We will explain the syntax later
*   DESCRIPTION -
    A fuller description based on the SYNOPSIS including a description of all the options to the command.
    This section may also include some example usage or extra details of how the command works.
*   AUTHOR -
    Who wrote the code to implement the command
*   REPORTING BUGS -
    How to report bugs
*   COPYRIGHT -
    Copyright information 
*   SEE ALSO -
    Points out other commands that you might find useful or which complement the current command.
    Also points to alternative sources of information.

#### How to Read the Synopsis

Here is the is synopsis for the 'cp' command:

~~~
SYNOPSIS
   cp [OPTION]... [-T] SOURCE DEST
   cp [OPTION]... SOURCE... DIRECTORY
   cp [OPTION]... -t DIRECTORY SOURCE...
~~~

This tells the reader that there are 3 ways to use the command. We will dissect
the first one:

~~~
cp [OPTION]... [-T] SOURCE DEST
~~~

This means the command `cp` followed by `[OPTION]...`. When something is in
square brackets, it means it can be left out. When something is followed by
`...` it means it can be repeated. Thus, `[OPTION]...` means any number of
options, including zero. `[-T]` means you can use the options `-T` here if you
wish. Then you add a SOURCE (the source file or directory) and a DEST (the 
destination file or directory).

The variants can be read in similar ways. Note that to use the last one, the
`-t` option is mandatory. In fact, this option switches the order of arguments
so that the destination immediately follows it.

#### Finding Help on Specific Options

The man output can serve as a quick reference for an option by searching. For
example, we may want to know more about the `-t` option mentioned above.

We can search the output of man with the slash key (/), followed by our search
query. To get information on `-t`, we would type `/-t` followed by RETURN.
Using 'n', we can navigate to the _next_ match until we find the detailed
information we need:

~~~
-t, --target-directory=DIRECTORY
     copy all SOURCE arguments into DIRECTORY
~~~

This means that this option has the short form `-t` and the long form
`--target-directory` and that it takes an argument. Its meaning is to copy all
the SOURCE arguments into DIRECTORY. Thus, we can give the destination
explicitly instead of relying on having to place the directory at the end.

#### What If I Cannot Remember the Command Name?

You can run man on itself

~~~
$ man man
~~~

and this will tell you that there is a '-k' option that does keyword searching
in the database of Linux commands on your current machine. So if I forget that 
the command is 'cp' but know it copies files I can use

~~~
$ man -k copy
~~~

to find command that relate to copying stuff. There is a quirky synonym for
this command, which might be simpler to remember (especially if you know
French):

~~~
$ apropos copy
~~~

This could produce a long list but you could filter it with the grep command as
follows

~~~
$ man -k copy | grep file
~~~

which helps a little and shows the 'cp' command at the top of the list.
