---
layout: lesson
root: ../..
title: Introduction
level: intermediate
---
Here's a typical working day for our robot.
The first thing it wants to do when it sits down in the morning is re-draw Figure 8 for the paper it's writing.
In order to do that, it needs to re-calculate the data (since it has some new numbers from a colleague back home).
In order to do *that*, it needs to recompile its statistics program, because it found and fixed a bug yesterday afternoon.
Oh, and it needs to update the version of Java installed on the machine it's using:
it has the right one on its laptop, but not on the machine it's using in the lab.
And it needs to re-install the graph-drawing tool that turns its data into a nicely-formatted chart suitable for publication.
That also requires a Java update, and it'll have to free up some disk space, too,
since *someone*'s MP3 collection is taking up 99.8% of available space.

And so it goes: everything it wants to do seems to require something else to be done first.
Until eventually, it finds himiself saying, "Right, so I need to shave the yak..."
We won't go that far with him, though&mdash;not today.

Instead, here's that set of jobs once again.
We can think of this set as a [graph](../../gloss.html#graph).

The things he needs to do&mdash;the [tasks](../../gloss.html#task)&mdash;are the graph's nodes.
The [dependencies](../../gloss.html#dependency) between the tasks are the graph's edges.
Since the robot can only do one thing at a time,
it needs to find an ordering on these tasks such that everything a given task depends on is done before it.

This pattern arises over and over again.
Each time we collect new data, we need to recalculate our summary statistics.
Each time our source files change, we need to recompile our program
(if, that is, we're using a compiled language like Java, C++, or Fortran).
And when someone writes a new paper, or receives an award, we need to update our research group's web site.
If there are more than a dozen or so tasks, it can be hard or impossible to manually keep track of what depends on what,
and what is up-to-date and what isn't, i.e., what's been done and what still needs to be done.

This is where tools like Make come in.
One of the fundametal rules of computing is that anything worth repeating is worth automating.
If we need to do the same tasks over and over again,
we should use a [build manager](../../gloss.html#build-manager) to handle the details.

We describe dependencies in a [build file](../../gloss.html#build-file),
which is usually just a plain text file in some specialized format.
We also describe how to update things,
i.e., what commands to run when something's dependencies have been satisfied and it's ready to be refreshed itself.
And that's all: the build manager handles everything else.
In particular, it keeps track of what's up to date, and what's ready to be updated.

The most widely used build manager on Unix and its derivatives is called Make.
And note that we said "most widely used", not "most popular".
Make was invented by a summer intern at Bell Labs in 1975.
(He went on to become a vice president at IBM and Google, which shows you how far a good program can take you.)
Over 35 years, Make has grown into a little programming language.
A very cryptic little language, without a debugger, whose conventions and rules only make sense if you understand the Unix shell.

The good news is, GNU Make (the de facto standard version of Make) is fast, free, and well-documented.
And many other tools know how to work with Make.
In particular, many integrated development environments can manage Make's build files more or less automatically,
shielding users from the ugly details.

In this chapter, we'll look at Make's basic features, and a few of its advanced facilities as well.
A companion lecture to this one explores a newer build manager called SCons.
It is more powerful and more flexible than Make, but isn't nearly as widely used (yet).
Java users should also look at Apache Ant, the standard build manager for Java.
It hides many of the platform-specific details that bedevil Make, but requires users to write XML files to get things done.
