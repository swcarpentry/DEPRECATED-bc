---
layout: lesson
root: ..
title: Motivation and Demotivation
---
## Relevance and Empowerment

Discussion of the practical implications of learning concepts
brings us to our next big idea:
people learn best when they care about the topic *and* believe they can master it.
Neither fact is particularly surprising,
but their practical implications have a lot of impact on what we teach
and the order in which we teach it.

First, most scientists don't actually want to program.
They want to do scientific research,
and programming is a tax they have to pay along the way.
They don't care how hash tables work, or even that hash tables exist;
they just want to know how to process data faster.
We therefore have to make sure that everything we teach is useful right away,
and conversely that we don't teach anything just because it's "fundamental".

Second,
believing that something will be hard to learn is a self-fulfilling prophecy.
This is why it's important not to say that something is easy:
if someone who has been told that tries it,
and it doesn't work,
they are more likely to become discouraged.

It's also why installing and configuring software is a much bigger problem for us
than experienced programmers like to acknowledge.
It isn't just the time we lose at the start of boot camps
as we try to get a Unix shell working on Windows
or set up a version control client on some idiosyncratic Linux distribution.
It isn't even the unfairness of asking students to debug things
that depend on precisely the knowledge they have come to learn,
but which they don't yet have.
The real problem is that every such failure reinforces the belief that computing is hard,
and that they'd have a better chance of making next Thursday's conference submission deadline
if they kept doing things the way they always have.

For these reasons,
we have adopted a "teach most immediately useful first" approach.
Imagine a grid whose axes are labelled "mean time to master" and "usefulness once mastered".
Everything that's quick to master, and immediately useful should be taught first;
things in the opposite corner
that are hard to learn and have little near-term application
don't belong in this course.

Many of the foundational concepts of computer science,
such as computability,
inhabit the "useful but hard to learn" corner of the grid described above.
This doesn't mean that they aren't worth learning,
but if our aim is to convince people that they *can* learn this stuff,
and that doing so will help them do more science faster,
they are less compelling than things like automating repetitive tasks.

And note:
any useful estimate of how long something takes to master must take into account
how frequent failures are
and how much time is lost to them.
For example,
editing a text file seems like a simple task,
but most GUI editors save things to the user's desktop or home directory.
If people need to run shell commands on the files they've edited,
a substantial fraction won't be able to navigate to the right directory without help.

## Mindset, Stereotype Threat, and Diversity

Carol Dweck and others have found that if you tell people that ability at some task is intrinsic
(i.e., that you either have the gene or you don't),
*everyone* does worse, including the supposedly advantaged.
The reason is that if they don't get it at first,
they figure they just don't have the gene,
which biases future performance.
On the other hand,
if people believe that a skill will improve with practice,
everybody does better on average.

A related concept is stereotype threat.
In brief,
reminding people of negative stereotypes,
even in subtle ways,
increases their nervousness and therefore their likelihood of failure.

The biggest negative stereotypes in computing are gender-related.
Depending on whose numbers you trust,
only 12-18% of programmers are women,
and those figures have actually been getting worse over the last 20 years.
There are many reasons for this
(see *Unlocking the Clubhouse* or the Whitecraft and Williams paper in the reading for discussion),
but as far as Software Carpentry goes,
the most important thing is to avoid both overt and covert reminders of social stereotypes of programmers.

We try to act on this in two ways.
First, we repeatedly emphasize that practice makes perfect.
We also code live in front of our learners instead of using slides,
so that they can see us make mistakes.
Doing this gives them permission to make mistakes too:
after all, if we're not perfect, they can't be expected to be either.
Having to type things in ourselves also slows us down,
so that learners are less likely to fall behind.

Second,
we have found that learners get more out of boot camps if they attend in groups,
e.g.,
if entire lab groups come,
or if attendees are drawn from the same (or closely-related) disciplines.
Doing this ensures that everyone in the room knows in advance
that they're going to be with at least a few people they trust,
and that they aren't going to be the only one in the room who doesn't get it.
It also helps after the workshop:
if someone comes on their own and then returns to their lab,
they have to roll a big rock up a steep hill all by themselves.
If they come with their labmates,
on the other hand,
they can work together after the bootcamp to implement what they've learned.

## Media-First Computation

FIXME: Guzdial and Ericson's work
