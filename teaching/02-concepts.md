---
layout: lesson
root: ..
title: Cognitive Development and Concept Mapping
---
## What Does It Mean to Understand Computing?

Our next problem is to figure out
what we want people to understand when we're finished teaching.
"How to write a loop" and "how to fetch records from a database"
are the practical skills we want them to be able to apply,
but these aren't useful on their own.
Without some higher-level understanding of computing,
people will still be stuck in a tweak-and-pray world.

Most definitions of what it means to understanding computing aren't very useful.
In particular,
the term "computational thinking" has been adopted by so many people that it now means little more than,
"Whatever the speaker thinks is important about computing."
For a better answer,
we need to turn to Mark Guzdial,
who has been studying computing and education for almost two decades.

[His answer](http://computinged.wordpress.com/2012/05/24/defining-what-does-it-mean-to-understand-computing/)
depends on two definitions.
First,
a *mental model* is a person's internal mental representation of something in the real world.
My mental model of how airplanes fly,
for example,
includes things like lift, thrust, banking, and fuel consumption.
It isn't physically accurate,
but it's close enough that I can predict what planes can and can't do well enough for my everyday needs.

Second,
a *notional machine* is the general properties of the computer that one is learning to control.
For example,
a notional machine for basic Python includes the idea of variables as sticky notes attached to values,
while one for C includes the notion of variables as names for locations in memory.
Neither is accurate,
but both are useful.

Given these definitions,
we can say that,
"To understand computing is to have a robust mental model of a notional machine."
In other words,
someone understands computing when their mental model of what the computer is doing
allows them to (more or less) predict how computers will behave.
For example,
if someone understands why this bit of code:

~~~
def change():
  var = 2

var = 1
change()
print var
~~~

prints 1 instead of 2,
then they understand call stacks in Python well enough to do most things,
even if they can't explain what a [closure](https://en.wikipedia.org/wiki/Closure_%28computer_programming%29) is
or define the term "lexical scoping".

Coming at it from another direction,
if someone has a robust mental model of a notional machine,
they can debug problems systematically rather than simply making random changes until it appears to work.
