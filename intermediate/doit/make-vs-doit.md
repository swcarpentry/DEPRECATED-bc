---
layout: lesson
root: ../..
title: A comparison of Make and doit
level: intermediate
---

## Why learn doit and not Make?

Software Carpentry has lessons covering two different build
tools: Make and doit. Which should you choose when starting
to automate your data analysis pipelines? This page is intended
to give a brief overview of the benefits and drawbacks of each
tool.

### 1. Popularity

Make is by far the most popular build tool currently in use. This
is probably it's greatest advantage. Learning how to write 
Makefiles will help you to reuse some pipelines written by other
people. It will also be a big help if and when you need to
install software from source code in a compiled language. If you
have ever installed something by typing some combination of
`./configure; make; make install` then you have already been
using Makefiles without knowing. Doit, on the other hand,
has a much smaller user base. If you learn doit, you may have
more difficulty finding examples, and less people will be able
to help you on forums or sites like Stack Overflow.

### 2. Availability

Make is available by default on most unix machines (i.e. linux 
and mac). Doit, on the other hand, is highly unlikely to be
pre-installed on any machine you might use. Doit is fairly
easy to install if you have root privileges, but suffers from 
all the problems of python packaging if you do not. On the other
hand, if you are a windows user, neither doit nor Make will
be available by default and doit may well be considerably easier
to install.

### 3. Readability

Doit is much more verbose and readable than Make. Make has a lot
of non-obvious automatic variables such as $^ and $@, whereas
automatic variables in doit are probably clearer, the equivalents
being %(dependencies)s and %(targets)s. Any if/then control flow
or parameters that you may wish to add are also likely to be
much more readable in doit/python than in make - especially if
the person reading does not know either language.

### 4. Self-Documentation

Doit implements a `list` command, which will list all the
available targets on the command line. If you write docstrings
for your tasks, doit will also expose these as help pages for
each target. By default, neither Make nor doit does a good job
of listing possible parameters. However, with a little work
doit can be combined with other python libraries such as
argparse, which allows for a much more usable and self-
documenting user interface

### 5. Dry Runs

Make can be run in "dry run" mode, where it won't actually do
anything, but will instead simply print out the commands
that would be executed. Doit does not have an equivalent of
this command.

### 6. Python

If you already know python, learning doit ought to be much
quicker than learning Make. Additionally, if you work mainly
in python, automating your pipelines with doit would reduce
the cognitive overload associated with switching between
programming languages.

### 7. Debugging

Doit can make use of any of the available python debuggers.
You can also run doit with the --pdb option to automatically
drop into a pdb shell when any errors occur. There are 
debuggers available for Make, but they are not especially 
easy to use.

### 8. Extensibility

Doit is very lightweight, and very extensible. You can easily
add your own commands, your own ways of determining whether
tasks should be re-run, your own ways to report the results
of a pipeline run, and much more.

### 9. The doit Database

Make decides whether a file has changed if the timestamp
on a dependency is newer than the timestamp on a target.
Doit decides a file has changed if it detects a change
in the files MD5 hash (i.e. only if the actual contents
of the file changes). This eliminates un-needed re-
calculation of some target files. On the other hand, in
order to do this doit needs to store file MD5 hashes in
a database file. If you accidentally delete this file,
or if you run doit in such a way that it can no longer
find the database file, it will be unable to tell which files 
have changed, and will default to re-running all of your 
tasks.

