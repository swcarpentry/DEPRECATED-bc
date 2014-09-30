---
layout: lesson
root: ..
title: Pedagogical Content Knowledge and Cognitive Load
---
Another way of looking at what teachers (ought to) know
divides knowledge into three parts:

- *content knowledge*, such as the "what" of programming;
- *general pedagogical knowledge*, i.e., an understanding of the
  psychology of learning; and
- the *pedagogical content knowledge* (PCK) that connects the two,
  which is things like when to teach the notion of a call stack
  and what examples to use when doing so.

This training course will focus on general pedagogical knowledge,
and will assume you know as much as you need to about basic programming (our content knowledge).
We *wish* the course could include a lot more material on PCK than it does,
but computing education lags 20-30 years behind fields like mathematics
when it comes to assembling and organizing that.
What Software Carpentry has is in [our teaching guide]({{page.root}}/novice/teaching/index.html);
additions and corrections would be very welcome.

## Cognitive Load

Another psychological theory that we use in our teaching is that of cognitive load.
The theory behind it is sometimes criticized,
but instruction based on it has been proven effective,
and it's a good framework for tying together several other ideas about learning.

In brief,
cognitive load theory holds that
people's brains are dealing with three kinds of load when they're learning:

- *Intrinsic* load is what they have to keep in mind
  in order to carry out a learning task.
- *Germane* load is the (desirable) mental effort required
  to create linkages between new information and old
  (which is one of the things that distinguishes learning from memorization).
- *Extraneous* load is everything else that distracts or gets in the way.

The key idea is pretty simple:
eliminating extraneous cognitive load accelerates learning.
The hard part is to figure out what's extraneous
(which is why the theory is often called unfalsifiable),
but research over the last three decades has identified a few factors.
One example is the work by Richard Mayer and others on the split-attention effect.
Correlating linguistic, auditory, and visual streams of information takes cognitive effort:
the brain can't help but check that it's getting the same information from those channels.
Learning is therefore more effective
when the same information is *not* being presented simultaneously in two different channels.
For example,
audio narration with on-screen captions is harder to learn from than either on its own;
speech and images is more effective *without* captions,
at least for native speakers of the language.

Second,
searching for a solution strategy is itself a large cognitive load.
This load can be reduced by giving learners worked examples,
i.e.,
by showing them a problem and a detailed step-by-step solution.
To maximize their impact,
worked examples should immediately be followed by a series of faded examples:
exercises in which learners are presented with a problem and a solution
in which some parts are left blank for them to fill in.
For example,
after having had this code explained to them:

~~~
def get_word_lengths(words):
    word_lengths = []
    for item in words:
        word_lengths.append(len(item))
    return word_lengths

print word_lengths(['hello', 'world'])
~~~

learners could be asked to fill in the blanks in:

~~~
def word_lengths(words):
    word_lengths = ____
    for item in range(len(____)):
        word_lengths.append(len(____))
    return word_lengths
~~~

After being shown:

~~~
doubles = [2 * x for x in items]
~~~

learners could be asked to get a list of lengths:

~~~
lengths = [____ for ____ in words]
~~~

and then put the ideas together in this template:

~~~
def word_lengths(____):
    return [____]
~~~

Faded examples are less intimidating than a blank screen.
In particular,
learners are much less likely to feel that they don't even know where to start.
They also encourage learners to think about the similarities and differences between various approaches,
which helps shape the conceptual categories we want them to form.
