Broadly speaking, teachers need three kinds of knowledge:

- *content knowledge*, such as the "what" of programming;
- *general pedagogical knowledge*, i.e., an understanding of the
  psychology of learning; and
- the *pedagogical content knowledge* that connects the two, which is
  things like when to teach the notion of a call stack and what
  examples to use when doing so.

In this training course,
we will assume you know enough about programming to teach it,
and will talk about general pedagogical principles.
We hope that this will help you make sense of the pedagogical content knowledge
that is implicitly presented in everything that has come before.

## What Does It Mean to Understand Computing?

Our starting point is the most basic question of all:
what are we trying to accomplish?
What do we want people to understand when we're finished teaching?
"How to write a loop" and "how to fetch records from a database"
are the practical skills we want them to be able to apply,
but these aren't useful on their own.
Without some higher-level understanding of what computing *is*,
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
a mental model is a person's internal mental representation of something in the real world.
My mental model of how airplanes fly,
for example,
includes things like lift, thrust, banking, and fuel consumption.
It isn't physically accurate,
but it's close enough that I can predict what planes can and can't do well enough for my everyday needs.

Second,
a notional machine is the general properties of the computer that one is learning to control.
For example,
a notional machine for basic Python includes the idea of variables as sticky notes attached to values,
while one for C includes the notion of variables as names for locations in memory.
Neither is completely accurate,
but both are workable.

Given these definitions,
we can say that,
"To understand computing is to have a robust mental model of a notional machine."
In other words,
someone understands computing when their mental model of what the computer is doing
allows them to (more or less) predict how computers will behave.
For example,
people don't need to know how pipes work in the shell in order to understand the shell;
what they need to know is how pipes behave.
One way to test this is to check that they understand the difference between:

~~~
$ head -10 data.txt | wc
~~~

and:

~~~
$ head -10 | wc data.txt
~~~

If they can explain why the answers are different,
and why the second will never finish if left undisturbed,
then they understand pipes.
Putting it another way,
if someone can identify the mistakes in something,
and debug something systematically rather than simply making random changes until it appears to work,
they probably understand it.

## There Are No Blank Slates

If understanding is where we're going, what's our starting point?
The answer is not "a blank slate" because there is no such thing.
Everyone has some mental model of the world, even at birth:
our visual cortex is wired for edge detection,
and other parts of our brain are wired for face recognition,
so we cannot help but put things in those mental boxes.

This fact explains the key difference between novice and competent practitioners in any field.
To paraphrase Tolstoy,
everyone who understands something understands it the same way,
but everyone who misunderstands,
misunderstands differently.
The reason is that by definition,
novices don't yet have the right conceptual categories for what they're learning,
so they are putting new knowledge into old boxes,
and those boxes are the wrong ones for this domain.
Everyone who has taught novices has seen first-hand that their questions often don't make sense;
this is the reason.

As an example of this idea's practical implications,
our introduction to the Unix shell only introduces 15 commands in two and a half hours.
Ten minutes per command may seem glacially slow,
but that's because our real goal is
to help learners construct the mental model and notional machine that these commands fit into.
That model includes things like:

- If you repeat things manually you'll eventually get them wrong,
  so let the computer repeat things for you by using tab completion and the `history` command.
- Lots of little tools, freely combined,
  are more productive than a handful of "kitchen sink" programs.

These two example illustrate something else as well.
Learning consists of more than "just" conveying facts and techniques:
creating linkages between concepts is as least as important.
Telling people that they shouldn't repeat things,
and that they should try to think in terms of little pieces loosely joined,
both set the stage for the introduction of functions in our programming lesson;
explicitly referring back to the shell lesson at that point helps solidify both.

### If You Use Robots To Teach, You Teach People To Be Robots

The transition from novice to competent practitioner
is no more and no less than the construction of correct (enough) categories,
i.e.,
the construction of a new mental model of this new intellectual domain.
The goal of education for novices is therefore
to help them form the right categories.
Until they've done that,
trying to impart "mere information" just confuses them.

This is one of the reasons most online courses don't work well.
Replaying recorded content is as ineffective for most learners as watching TV,
or as a professor standing in front of 400 people in a lecture hall,
because neither can intervene to clear up specific learners' misconceptions.
Some people happen to already have the right conceptual categories for a subject,
or happen to form them correctly early on;
these are the ones who stick with most massive online courses,
but many discussions of the effectiveness of such courses ignore this survivor bias.

### The Documentation Is Useless

This is also one of the reasons software documentation is so often frustrating.
Reference material is opaque to someone who doesn't know what they're looking for,
such as a novice who doesn't yet have a mental map of the domain.
On the other hand,
tutorials meant to help people build such a map
are too slow and too diffuse for people who already have one.
It is possible to craft something that serves both communities,
but it's often simpler to address their needs separately.

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
Imagine a grid whose axes are labeled "mean time to master" and "usefulness once mastered".
Everything that's quick to master, and immediately useful should be taught first;
things in the opposite corner
that are hard to learn and have little near-term application
don't belong in this course.

Many of the foundational concepts of computer science,
such as computability,
inhabit the "useful but hard to larn" corner of the grid described above.
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
they have to roll a big rock up a steep hill all by themself.
If they come with their labmates,
on the other hand,
they can work together after the bootcamp to implement what they've learned.

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

and then work up to:

~~~
def word_lengths(words):
    return [____ for ____ in ____]
~~~

Faded examples are less intimidating than a blank screen.
In particular,
learners are much less likely to feel that they don't even know where to start.
They also encourage learners to think about the similarities and differences between various approaches,
which helps shape the conceptual categories we want them to form.

Scaling Instruction
-------------------

We said earlier that clearing up novices' misconceptions is more important than imparting "mere information".
However,
we also said that when novices misunderstand something,
they can all misunderstand it differently.
How can one instructor clear up the different misconceptions of forty or four hundred people at a time?

The best solution we have so far is called peer instruction.
Originally developed by Eric Mazur,
it has been studied extensively in a wide variety of contexts,
including programming.
When it is used,
the basic learning cycle is typically something like this:

1. The instructor presents a multiple-choice question.
2. Learners choose an answer
   (typically using clickers or by holding up multi-colored cards).
3. Learners discuss their answers in small groups (typically 3-4 people) and re-vote.
4. The instructor presents and explains the correct answer.
5. Learners re-group to discuss the solution.

The second discussion stage is crucial,
because that's where group members clear up each other's misconceptions.
This part of the cycle is what makes the technique scale:
with three or four people per group,
the odds are good that at least one person will understand what a partner has misunderstood,
and be able to clarify it for them.

We don't use peer instruction in our workshops,
partly because it depends on people having looked at material before they come to class
(which most won't actually do),
and partly because learning how to take part in peer instruction is extraneous cognitive load
for people who are also trying to learn how to program.
What we've done instead is require learners to do exercises in pairs,
i.e.,
have two (or three) people share a laptop and talk to each other as they do each exercise.
This practice is called pair programming,
and many studies have shown its effectiveness in industry.
It works because weak learners get individualized explanations,
while stronger ones learn from explaining things.
It's particularly effective if pairs are re-formed periodically,
i.e.,
if people are required to sit beside someone new after every break
so that weak-weak and strong-strong pairs don't persist.
(On the other hand,
this can make some people more inhibited,
since eventually they'll be paired with a stranger.)

## Summing Up

Teaching is a performing art,
but anyone who wants to do it well can draw on an enormous body of research.
The most important thing is to realize that while practice makes perfect,
feedback makes perfect faster:
teaching with a partner is the best way to improve what you do and how you do it.
