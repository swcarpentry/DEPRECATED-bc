---
layout: lesson
root: ../..
title: Automating Tasks with Make
level: intermediate
---

Originally invented to manage compilation of programs, `make` can be
used to automatically execute sequences of commands (programs) and
update any set of files that depend on another set of files. This
makes it a good solution for many data analysis and data management
problems, including the generation of images from data.

Programs that `make` will execute are described by a text file (almost
always called `Makefile`) containing a list of commands, called
*rules*, and the files that they create, called *targets*. The rules
describe how files depend on each other, and how to update out-of-date
files. They can be very specific, useful to generate one file, or
general, specifying a pattern for creating a certain class of files.

Topics
------
*   [Introduction](00-intro.html)
*   [Basic Tasks](01-basics.html)
*   [Automatic Variables and Wildcards](02-automatic-variables.html)
*   [Patterns](03-patterns.html)
*   [Macros](04-macros.html)

See Also
--------
*   [Instructor's Guide](guide.html)
*   [Reference](ref.html)

Resources
---------
*   \*.mk: Makefiles used as examples in notes.
