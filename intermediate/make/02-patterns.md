---
layout: lesson
root: ../..
title: Patterns
level: intermediate
---
Let's go back to our paper and look at another part of our dependency graph.
`summary-1.dat` depends on all of the files `data-1-1.dat`, `data-1-2.dat`, and so on.
The number of files isn't fixed: there could be one, a dozen, or a thousand.
Writing a rule for exactly three files is easy&mdash;we just have one target and multiple prerequisites on a single line:

    # multiple.mk

    summary-1.dat : data-1-1.dat data-1-2.dat data-1-3.dat
            stats.py summary-1.dat data-1-1.dat data-1-2.dat data-1-3.dat

But how do we generalize that to any number of files?
And how can we get rid of the repeated filenames?
Writing `data-1-1.dat data-1-2.dat data-1-3.dat` twice is just asking for trouble:
sooner or later, we'll add a file to one line but forget to update the other.
We need a way to express the idea, "All the files named `data-1-something.dat`,"
even when we don't know in advance how many of these there will be.
We'd also like to figure out what to do about `figure-2.svg` and the files it depends on.
The rules are exactly the same as those for `figure-1.svg` and its prerequisites;
duplicating them is just asking for trouble.

Let's start with the case of three files `data-1-1.dat`, `data-1-2.dat`, and `data-1-3.dat`.
As we have seen,
it's easy to write a Make rule to update `summary-1.dat` whenever any of these or the `stats.py` script change.

We'd like to do better, though, so let's replace the action in the rule:

    # target-variable.mk

    summary-1.dat : data-1-1.dat data-1-2.dat data-1-3.dat
            stats.py $@ data-1-1.dat data-1-2.dat data-1-3.dat

Instead of naming `summary-1.dat` in the rule's action, we use the rather cryptic shorthand `$@`.
This is one of Make's [automatic variables](../../gloss.html#automatic-variable),
and it means "the target of the current rule".
In this rule, for example, it means `summary-1.dat`.
(And no, there isn't a more readable long form of the name: it's just another of Make's many warts.)

Using `$@` instead of repeating the target's name shortens our rule somewhat,
but writing the many prerequisite filenames twice is still redundant.
Let's fix that by replacing our shortened rule command like this:

    # variables.mk

    summary-1.dat : data-1-1.dat data-1-2.dat data-1-3.dat
            stats.py $@ $^

`$^` is another automatic variable: it means "all the prerequisites of this rule".
In this case it's the three raw data files,
so when Make expands the variables in `stats.py $@ $^`,
we get back our original command.

There are other automatic variables as well:
for example, `$<` means "the first prerequisite in the list",
and `$?` means "all prerequisites that are out of date".
Don't worry if you can't remember them:
everyone except the most passionate Make user writes them on a sticky note and puts it on their terminal.

Using the automatic variables `$@` and `$^` eliminates the redundancy in our rule,
but doesn't solve the problem of handling an arbitrary number of prerequisite filenames.
We expect to have more than three data files before this project is done, and as we said before,
we don't want to have to rewrite our Makefile each time we run our experiment.
What we really want is something like the shell's `*` wildcard, which matches any number of characters:

    # wildcard.mk

    summary-1.dat : data-1-*.dat
            stats.py $@ $^

This actually works:
if use `data-1-*.dat` as the rule's prerequisite, it behaves just like the corresponding shell wildcard.
When we do this, we *must* use `$^` to refer to the rule's prerequisites in the action:
we don't know exactly what filenames will match,
so we have to rely on Make to put them in an automatic variable for us on a rule-by-rule basis.

Here are our dependency tree and our entire Makefile so far:

    paper.pdf : paper.wdp figure-1.svg figure-2.svg
            wdp2pdf $<

    figure-1.svg : summary-1.dat
            sgr -N -r $@ $^

    figure-2.svg : summary-2.dat
            sgr -N -r $@ $^

    summary-1.dat : data-1-*.dat
            stats.py $@ $^

    summary-2.dat : data-2-*.dat
            stats.py $@ $^

There is still some redundancy:
we have exactly the same logical rules for our two data series,
but have to write them down separately because the '1' and '2' in their names are different.

We'll see how to fix this in the next section.
Before then, though, we have one more problem to address.
Our existing Makefile doesn't capture the fact that `summary-1.dat` and `summary-2.dat`
depend on `stats.py` as well as on their corresponding raw data files.
We could try to fix this by adding `stats.py` to their prerequisite lists:

    paper.pdf : paper.wdp figure-1.svg figure-2.svg
            wdp2pdf $<

    figure-1.svg : summary-1.dat
            sgr -N -r $@ $^

    figure-2.svg : summary-2.dat
            sgr -N -r $@ $^

    summary-1.dat : stats.py data-1-*.dat
            stats.py $@ $^

    summary-2.dat : stats.py data-2-*.dat
            stats.py $@ $^

If we do this, though, `stats.py` will appear in the value of the automatic variable `$^` for those two rules.
This means that when we run `stats.py`,
our command line will be `stats.py summary-1.dat stats.py data-1-1.dat data-1-2.dat` and so on,
i.e., we'll be telling `stats.py` to process itself as a data file, which is almost certainly a bad idea.
We could "fix" this by having `stats.py` ignore files that end in `.py`, but it would be an ugly hack.

A second option would be to move the dependency down, and pretend that the raw data files depend on `stats.py`:

    figure-2.svg : summary-2.dat
            sgr -N -r $@ $^

    summary-1.dat : data-1-*.dat
            stats.py $@ $^

    summary-2.dat : data-2-*.dat
            stats.py $@ $^

    data-1-1.dat : stats.py
            touch $@

    data-1-2.dat : stats.py
            touch $@

This is called a [false dependency](../../gloss.html#false-dependency).
The raw data files don't really have to be updated when `stats.py` is changed,
but with this false dependency in our Makefile,
Make will update the timestamps on the raw data files when `stats.py` changes,
which will in turn trigger an update of the summary files.

False dependencies do solve some problems, but not this one:
if we go down this road, we have to list all our raw data files explicitly once again, which is what we're trying to avoid.
Our third option is
to add additional rules for `summary-1.dat` and `summary-2.dat`
that add `stats.py` as a prerequisite,
but don't have any actions:

    paper.pdf : paper.wdp figure-1.svg figure-2.svg
            wdp2pdf $<

    figure-1.svg : summary-1.dat
            sgr -N -r $@ $^

    figure-2.svg : summary-2.dat
            sgr -N -r $@ $^

    summary-1.dat : data-1-*.dat
            stats.py $@ $^

    summary-2.dat : data-2-*.dat
            stats.py $@ $^

    summary-1.dat : stats.py
    summary-2.dat : stats.py

When Make sees multiple rules for the same target,
it uses the union of all those rules' prerequisites as the target's actual set of prerequisites.
However, the automatic variable `$^` in the rule is still just that rule's prerequisite list.
It's a bit of a hack, but it means that our command line has exactly what we want it to have.
