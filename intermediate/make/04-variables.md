---
layout: lesson
root: ../..
title: Variables
level: intermediate
---
This chapter describes yet another way to make Makefiles rules more
generic. Make *variables* can contain be used to specify any part of a
rule (the target, the prerequisites, the commands to run, parts of
file names, etc.). They can also be specified externally when invoking
`make`, adapting the rules without modifying the Makefile.

Objectives
----------
* show the general syntax for using Makefile variables
* show how variables can be defined in the Makefile
* show how variables can be specified when invoking `make`

Just when we thought we were done writing our Makefile,
our supervisor reminded us that all papers must conform to the university's new style rules.
That means that `paper.pdf` has one more dependency: the official university style file `euphoric.wps`.
Unfortunately, on our laptop, that file lives in `C:\papers`,
but on the machine we use in the lab, it's in `/lib/styles`.

We could create a directory called `/lib/styles` on our laptop,
and put a copy of `euphoric.wps` there,
but the university also has a style guide for diagrams, which is in a file called `euphoric.fig`.
Once again, on our laptop, it's installed in `C:\papers`,
but it's in `/lib/styles` in the lab.
How should we handle this difference?

If we start with the Makefile we've written so far,
the brute-force approach is to just add the style files to our commands:

~~~
paper.pdf : paper.wdp figure-1.svg figure-2.svg
        wdp2pdf --style c:/papers/euphoric.wps $<
~~~
{:class="in"}

~~~
figure-%.svg : summary-%.dat
        sgr -N -r -s c:/papers/euphoric.fig $@ $^
~~~
{:class="in"}

~~~
summary-%.dat : data-%-*.dat
        stats.py $@ $^
~~~
{:class="in"}

~~~
data-*-*.dat : stats.py
        touch $@
~~~
{:class="in"}

There's some redundancy here, though: we are specifying the same directory twice.
And notice that we haven't explicitly listed `euphoric.wps` or `euphoric.fig`
as prerequisites of `paper.pdf`,
or of the two figure we're generating.
Some people would include them, just to be safe,
but it's more common not to list dependencies on "system" files
like libraries, style files, and so on.

Now, how do we handle the fact that these two paths need to be different when we're re-creating our paper in the lab?
The first option is to use copy and paste, and write two completely separate Makefiles.
What we really mean, though, is write *and maintain*, and that's why this is a bad idea.
As soon as we have two of anything, we'll eventually update one but forget to update the other.
Makefiles are already hard enough to debug;
any "solution" that adds more complexity and risk isn't really a solution at all.

Our second option is to put everything in one Makefile, and then to comment out the bits intended for the machine we *aren't* on,
but this is also problematic.
First, we have to make sure we always comment and uncomment lines consistently.
If we uncomment the line for creating the paper on our laptop, for example, but forget to uncomment the line for building the figures,
we're going to have another debugging headache.

Commenting and uncommenting lines also makes life more difficult for our version control system.
If we update our Makefile from version control, then change the commenting on a few lines,
the version control system will want to save those changes in the repository the next time we commit.
We probably don't actually want to do that, since it would mean that the next time we updated on the other machine,
its Makefile would be overwritten.

The third option&mdash;the right one&mdash;is to refactor our Makefile to make the problem go away entirely.
We can do this by defining a [variable](../../gloss.html#variable), just as we would define a constant or variable in a program.
Here's our Makefile with a variable defined and used:

~~~
# with-macro.mk
STYLE_DIR=c:/papers/

paper.pdf : paper.wdp figure-1.svg figure-2.svg
        wdp2pdf --style ${STYLE_DIR}/euphoric.wps $<

figure-%.svg : summary-%.dat
        sgr -N -r -s ${STYLE_DIR}/euphoric.fig $@ $^

summary-%.dat : data-%-*.dat
        stats.py $@ $^

data-*-*.dat : stats.py
        touch $@
~~~
{:class="in"}

The definition looks like definitions in most programming languages:
the variable is called `STYLE_DIR`, and its value is `c:/papers/`.
To use the variable, we put a dollar sign in front of it (just as we would do in the shell) and wrap its name in curly or round brackets.
This tells Make to insert the variable's value, so that these two directory paths are what we want on our laptop.

This is certainly a step forward:
now, when we want to move our Makefile from one machine to another, we only have to change one definition in one place.
However, while we no longer have to worry about consistency,
we're still making changes to a file that's under version control that we *don't* want written back to the repository.

> #### Parenthesizing Variables in Make
>
> We have to put curly brackets or parentheses around a variable's name when we use it&mdash;we can't just write `$VARIABLE`.
> If we do, `make` will interpret it as `$V` (a reference to the variable `V`) followed by "ARIABLE".
> Since we probably don't have a variable called `V`, `$V` will expand to the empty string,
> so `$VARIABLE` without parentheses will just be "ARIABLE".
> Why?
> To make a long story short, it's another wart left over from history.
> Almost everyone trips over it occasionally, and as with other bugs, it can be very hard to track down.

Using variables for program parameters
--------------------------------------

It's common practice to use variables to define all the flags that tools need,
so that if a tool is invoked in two or more actions,
it's passed a consistent set of flags.
Here, for example, we're defining `STYLE_DIR` to point to the directory holding our style files,
then using that definition in two other variables:

~~~
# with-lots-of-macros.mk

STYLE_DIR=c:/papers/
WDP2PDF_FLAGS=--style ${STYLE_DIR}/euphoric.wps
SGR_FLAGS=-N -r -s ${STYLE_DIR}/euphoric.fig

paper.pdf : paper.wdp figure-1.svg figure-2.svg
        wdp2pdf ${WDP2PDF_FLAGS} $<

figure-%.svg : summary-%.dat
        sgr ${SGR_FLAGS} $@ $^

summary-%.dat : data-%-*.dat
        stats.py $@ $^

data-*-*.dat : stats.py
        touch $@
~~~
{:class="in"}

The first, `WPD2PDF_FLAGS`,
is the single flag and argument we want to pass to the tool that turns our word processor file into a PDF.
The second, `SGR_FLAGS`, combines `STYLE_DIR` with a couple of other flags
to build the arguments for the tool that turns data files into SVG diagrams.

We are now ready to solve our original problem.
Let's move the definition of `STYLE_DIR`&mdash;the variable that changes from machine to machine&mdash;out of our main Makefile,
and into a Makefile of its own called `config.mk`:

~~~
# config.mk
~~~
{:class="in"}

~~~
STYLE_DIR=c:/papers/
~~~
{:class="in"}

We can then include that file in our main Makefile using Make's `include` command.
Our other variables and commands can then use the definition of `STYLE_DIR` just as if it had been defined in the main Makefile:

~~~
# with-include.mk
include config.mk
~~~
{:class="in"}

~~~
WDP2PDF_FLAGS=--style ${STYLE_DIR}/euphoric.wps
SGR_FLAGS=-N -r -s ${STYLE_DIR}/euphoric.fig
~~~
{:class="in"}

~~~
paper.pdf : paper.wdp figure-1.svg figure-2.svg
        wdp2pdf ${WDP2PDF_FLAGS} $<
~~~
{:class="in"}

~~~
figure-%.svg : summary-%.dat
        sgr ${SGR_FLAGS} $@ $^
~~~
{:class="in"}

~~~
summary-%.dat : data-%-*.dat
        stats.py $@ $^
~~~
{:class="in"}

~~~
data-*-*.dat : stats.py
        touch $@
~~~
{:class="in"}

Once we've tested this to make sure it works, we can copy `config.mk` to create two files that we'll put in version control.
The first, `config-home.mk`, defines `STYLE_DIR` for use on our laptop.
The second, `config-lab.mk`, defines it for use in the lab.
These two files are only changed when they need to be (i.e., when the style files move, or their names change).
We then copy one or the other on the machine we're using to create the file `config.mk`
that our main Makefile actually includes.

For example, here's what we have in the `paper` directory on our home machine when we do a fresh checkout from version control.
Along with our data files and the word processor file, we have our main Makefile and the two machine-specific configuration makefiles.
So we copy `config-home.mk` to create `config.mk`.
Meanwhile, when we check out in the lab, we copy `config-lab.mk` to create `config.mk`.
Our main Makefile is now happy in both cases because the file it's including now exists,
and has the right definition of `STYLE_DIR`.

We can also solve this problem by defining `STYLE_DIR` on the command line each time we run Make.
To do this, we set the variable on the command line when invoking `make`:

~~~
$ make STYLE_DIR=/lib/styles -f Makefile
~~~
{:class="in"}

This is often a bad idea, though.
We have to remember to type the definition each time,
and we have to type it *correctly* each time, or rely on shell history to retrieve the
right invocation.
The biggest problem is that it leaves no record in the Makefile itself of the settings used,
which makes life harder for other people who want to re-create our paper:
how do they know what to type?
