---
layout: lesson
root: ../..
title: Rules
level: intermediate
---
Our Makefile is better than it was, but still contains a lot of redundancy.
The rules for `figure-1.svg` and `figure-2.svg` are identical except for the '1' and '2' in their names,
as are the rules for `summary-1.dat` and `summary-2.dat`.

We'd like to "fold" the rules for the figures together for two reasons.
First, if we add a third figure, we don't want to have to duplicate rules a third time.
Second, if we ever want to change the way we generate figures, we'd like to make that change once, in one place:
if we have to make it in several places, the odds are good we'll forget one,
and then waste time trying to figure out why some of our commands aren't running.

The way to do this in Make is to use a [pattern rule](../../gloss.html#pattern-rule) to capture the common idea.
Here's our Makefile rewritten to use such a rule:

    # pattern-rule.mk

    figure-%.svg : summary-%.dat
            sgr -N -r $@ $^

    summary-1.dat : data-1-*.dat
            stats.py $@ $^

    summary-2.dat : data-2-*.dat
            stats.py $@ $^

    summary-1.dat : stats.py
    summary-2.dat : stats.py

In this rule, `%` is a wildcard.
When it is expanded, it has the same value on both sides of the rule:
if it matches '1' on the left, it must match '1' on the right as well.
`%` only means something to Make, though.
It doesn't have a value in the rule's action, which is handed off to the shell for execution.
So in the action, we have to use the automatic variables `$@` and `$^` as before.

Let's try running our modified Makefile:

    $ make -f pattern-rule.mk
    stats.py summary-1.dat data-1-1.dat data-1-2.dat data-1-3.dat

`summary-1.dat` is updated, but not `summary-2.dat` or either of the figure files.
The reason the other commands didn't run is that pattern rules don't create dependencies:
they just tell Make what to do *if* there's a dependency.
In other words, *if* Make decides it wants to create `figure-1.svg`, it can use our pattern rule,
but we still have to tell Make to care about `figure-1.svg`.
Let's do this by putting the rule for `paper.pdf` back in our Makefile:

    # use-pattern.mk

    paper.pdf : paper.wdp figure-1.svg figure-2.svg
            wdp2pdf $<

    figure-%.svg : summary-%.dat
            sgr -N -r $@ $^

    summary-1.dat : data-1-*.dat
            stats.py $@ $^

    summary-2.dat : data-2-*.dat
            stats.py $@ $^

    summary-1.dat : stats.py
    summary-2.dat : stats.py

Here, `paper.pdf` depends on `figure-1.svg` and `figure-2.svg`.
Make now knows that it needs these figures.
Since there aren't specific rules for them, it uses the pattern rule instead.

It's tempting to go one step further, and make `paper.pdf` depend on `figure-*.svg`:

    paper.pdf : paper.wdp figure-*.svg
            wdp2pdf $<

This doesn't work, though.
The reason is that the figure files may not exist when Make starts to run&mdash;after all, Make creates them.
In that case, `figure-*.svg` will expand to nothing,
so Make would mistakenly believe that `paper.pdf` depended only on `paper.wdp`.
This kind of bug can be very hard to figure out,
and while Make does have a debugger called [GMD](http://gmd.sourceforge.net/),
it's not an easy tool for beginners to use.

Our raw data files *do* always exist, though, so we can get rid of some more redundancy by folding these two rules into one
using the `*` wildcard:

    # all-patterns.mk

    paper.pdf : paper.wdp figure-1.svg figure-2.svg
            wdp2pdf $<

    figure-%.svg : summary-%.dat
            sgr -N -r $@ $^

    summary-%.dat : data-%-*.dat
            stats.py $@ $^

    summary-1.dat : stats.py
    summary-2.dat : stats.py

It's safe to do this because Make isn't responsible for creating `data-1-whatever.dat` and `data-2-whatever.dat`:
there's no possibility of the `*` missing things because it's evaluated when Make starts running.

Just as a reminder, the `%` is a Make wildcard:
it matches the same thing on the left and right side of a pattern rule.
`*` is a shell wildcard:
it matches zero or more characters in a filename when it's evaluated.

We cannot get rid of the last bit of redundancy by making `summary-%.dat` depend on `stats.py`.
Even with this pattern rule, the summary files only depend on the corresponding raw data files, not on `stats.py`.
The reason is that when Make sees two or more pattern rules that could match a filename,
it uses the first and ignores the other.
It's another wart, and another source of hard-to-find headaches in Makefiles.

If we really want to avoid making `summary-1.dat` and `summary-2.dat` depend on `stats.py` separately,
the only way is to go back to using false dependencies.
This Makefile tells Make to update the timestamps on the raw data files using `touch` whenever `stats.py` changes.
Doing this indirectly triggers the re-creation of the summary files&mdash;it does what we want, just in a roundabout way.

    # false-dependencies.mk

    paper.pdf : paper.wdp figure-1.svg figure-2.svg
            wdp2pdf $<

    figure-%.svg : summary-%.dat
            sgr -N -r $@ $^

    summary-%.dat : data-%-*.dat
            stats.py $@ $^

    data-*-*.dat : stats.py
            touch $@
