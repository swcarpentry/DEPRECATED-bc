---
layout: lesson
root: ../..
title: Basic Tasks
level: intermediate
---
To illustrate how Make works, here's the dependency tree for the paper that the robot is working on.
`paper.pdf` depends on `paper.wdp` (the raw word processor file),
and on `figure-1.svg` and `figure-2.svg`.
`figure-1.svg` depends on `summary-1.dat`,
which in turn depends on `data-1-1.dat`, `data-1-2.dat`, and so on,
while `figure-2.svg` depends on files with similar names.

In order to create `paper.pdf`, we have to run the command `wdp2pdf paper.wdp`.
For the purpose of this lecture, it doesn't matter what `wdp2pdf` actually does.
All we need to know is that if `paper.wdp` or either of the figure SVG's change, we need to re-run this command.

To create `figure-1.svg`, we run `sgr -N -r summary-1.dat` and send the output to `figure-1.svg`.
Again, it doesn't matter for now what the `sgr` command actually is.
What matters is that we need to run it whenever `figure-1.svg` is out of date,
i.e., whenever it is older than the `summary-1.dat` file it depends on.
Finally, in order to update `summary-1.dat`, we need to run our own little script, `stats.py`,
with all the files named `data-1-something.dat` as input.
We don't know in advance how many of these there will be: we could conceivably have dozens or hundreds of raw data files to summarize.

That little program `stats.py` adds one more wrinkle to our example problem.
We're constantly updating it as we think of new ways to process our raw data files.
We're also finding and fixing bugs more often than we'd like.
Each time it changes, we should probably update `summary-1.dat`,
just in case a new feature or bug fix changes the summary values.
We should therefore include `stats.py` in the list of things `summary-1.dat` depends on,
so that changes to `stats.py` will trigger recalculation of `summary-1.dat`.

This is all a bit much to digest at once, so let's look at the simplest piece.
How can we get Make to re-create `figure-1.svg` automatically whenever `summary-1.dat` changes?

Let's start by going into the directory containing the files we're using in the paper,
and use the `ls` command to get a listing of what's there.
The `-t` flag to `ls` tells it to list things by age, with the youngest file first and the oldest last:

    $ ls -t *.dat *.svg
    summary-1.dat    figure-1.svg

This listing tells us that our data file `summary-1.dat` is newer than the SVG file that depends on it,
so the SVG file needs to be re-created.
Using our favorite editor, let's create a file called `hello.mk` and put these three lines in it.

    # hello.mk
    figure-1.svg : summary-1.dat
            sgr -N -r summary-1.dat > figure-1.svg 

A configuration file for Make like this one is called a [Makefile](../../gloss.html#makefile).
The first line, starting with `#`, is a comment.
(Our comments should be more meaningful than just the name of the file.)
The second and third lines are a [rule](../../gloss.html#rule) that tell Make what we want to do.

The filename on the left of the colon in the first line is the [target](../../gloss.html#target) of the rule.
The rule tells Make how to update or re-create this file.
The target's [prerequisites](../../gloss.html#prerequisite)&mdash;the things it depends on&mdash;are listed to the right of the colon.
In our case, `figure-1.svg` only has one prerequisite, `summary-1.dat`.

The second line of the rule is its [action](../../gloss.html#action).
This tells Make what shell command or commands to run to bring the target up to date if it is older than any of its prerequisites.
This rule only has one command, but a rule can contain any number.

One thing to note is that the actions in rules *must* be indented with a single tab character.
Make will not accept spaces, or mixes of spaces and tabs.
(As we said in the introduction, it was written by a summer intern in 1975, and sometimes that shows.)

Now that we've created our Makefile, we can tell Make to obey its instructions by running `gmake` from the command line:

    $ gmake -f hello.mk
    sgr -N -r summary-1.dat > figure-1.svg

Many systems make `make` an alias for `gmake`,
so if the latter doesn't work for you, try the former name as well.
The arguments `-f hello.mk` tell Make that we want it to use the commands in the file `hello.mk`.
If we don't tell it what file to look in,
it looks for a file called `Makefile` in the current directory and uses that if it exists.

Make's output shows us that it has run the command we wanted it to.
It did this because at least one prerequisite for `figure-1.svg` was newer than `figure-1.svg` itself.
By default, Make uses the time a file was last modified as its age.
(Opening a file in an editor to view it doesn't change this timestamp, but any change to its contents will.)
Since `summary-1.dat`'s timestamp was younger than `figure-1.svg`'s,
Make ran the shell command we gave it and created a new version of `figure-1.svg`.

Let's run Make again:

    $ gmake -f hello.mk

This time, it doesn't execute any commands.
This happened&mdash;or didn't&mdash;because the target is newer than its prerequisites.
Since there's nothing to bring up to date, Make doesn't change anything.

If we were only allowed one rule per file,
Make wouldn't be any simpler than typing commands by hand or putting them in little shell scripts.
Luckily, Make allows us to put any number of rules in a single configuration file.
Here is another Makefile called `double.mk` with rules to re-create
both `figure-1.svg` and `figure-2.svg`.
These rules are identical except for the 1's and 2's in the filenames; we'll see later how to combine these rules into one.

    # double.mk
    figure-1.svg : summary-1.dat
            sgr -N -r summary-1.dat > figure-1.svg

    figure-2.svg : summary-2.dat
            sgr -N -r summary-2.dat > figure-2.svg

Let's pretend we've just updated our data files by running `touch *.dat`.
(The Unix `touch` command doesn't change the contents of files, but updates their timestamps as if they had been modified.)
Now, when we run Make, it re-creates `figure-1.svg` again&mdash;and then stops:

    $ touch *.dat
    $ gmake -f double.mk
    sgr -N -r summary-1.dat > figure-1.svg

Why wasn't `figure-2.svg` re-created?
The answer is that Make uses the first rule in the Makefile as its [default rule](../../gloss.html#default-rule).
Unless it's told otherwise, it only executes this rule.
If we want Make to rebuild `figure-2.svg`, we have to tell it so explicitly.
We use `-f double.mk` to tell Make what Makefile to use,
and then give it the name of the target we want it to handle:

    $ gmake -f double.mk figure-2.svg
    sgr -N -r summary-2.dat > figure-2.svg

Again, building things one at a time like this is slightly better than typing individual commands, but only slightly.
To get Make to build everything at once, we have to introduce a [phony target](../../gloss.html#phony-target).
This is just a target name that doesn't correspond to any actual file.
Since it doesn't actually exist, it can't ever be up to date, but other things can still depend on it.
Here's our third Makefile, `phony.mk`:

    # phony.mk

    all : figure-1.svg figure-2.svg

    figure-1.svg : summary-1.dat
            sgr -N -r summary-1.dat > figure-1.svg

    figure-2.svg : summary-2.dat
            sgr -N -r summary-2.dat > figure-2.svg

We've introduced a phony target called `all`, which depends on `figure-1.svg` and `figure-2.svg`.
Since there's no file called `all` in the current directory,
if we type `make all`,
Make will decide that the `all` target is out of date.
And since `all` depends on `figure-1.svg` and `figure-2.svg`,
Make will go and update them both, which is exactly what we want.

Let's `touch` our data files again, and run `make -f phony.mk all`.
Sure enough, Make runs the `sgr` command twice to re-create both figures:

    $ touch *.dat
    $ gmake -f phony.mk
    sgr -N -r summary-1.dat > figure-1.svg
    sgr -N -r summary-2.dat > figure-2.svg

One thing to note is that the order in which commands are executed is arbitrary.
Make could decide to update `figure-2.svg` first, rather than `figure-1.svg`,
because there's no dependency to respect between the two.
Make could also update them in parallel if it had more than one processor to use&mdash;we'll return to this idea later.

Something else this example shows us is that a single thing can be a target in one rule, and a prerequisite in others.
The dependencies between the files mentioned in the Makefile make up a directed graph.
In order for Make to run, this graph must not contain any cycles.
For example, if X depends on Y, Y depends on Z, and Z depends on X,
everything depends on something else, so there is nothing Make can build first.
If it detects a cycle in a Makefile, Make will print an error message and stop.
Unfortunately, whether or not a cycle exists depends on which files exist,
and Make's error message is usually not particularly informative.
