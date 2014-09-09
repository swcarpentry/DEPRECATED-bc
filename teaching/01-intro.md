---
layout: lesson
root: ..
title: Introduction
---
Let's start by defining a few terms.
The first is *educational psychology*,
which is the study of how people learn.
This touches on everything from the neuropsychology of perception and the mechanisms of memory
to the sociology of school systems
and the philosophical question of what we actually mean by "learning"
(which turns out to be pretty complicated once you start looking beyond
the standardized Western classroom).

What we know about educational psychology constrains teaching,
but doesn't dictate it.
In particular,
there are often several plausible approaches to teaching
that are consistent with what we know about how brains learn.
In order to decide between them,
we need to explore *instructional design*,
which is the study of how to create learning materials.
If ed psych is science,
ID is engineering:
it's how understanding is put into practice.

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

## Stages of Cognitive Development

The next piece of background we need is
how learners' mastery of a subject changes.
In the 1980s,
Patricia Brenner developed
[a five-stage model of skill acquisition](http://www.amazon.com/From-Novice-Expert-Excellence-Commemorative/dp/0130325228/)
that has been applied successfully in many fields.
For our purposes,
we can simplify that model to three stages:

- A *novice* is someone who doesn't know what she doesn't know.
  She doesn't yet have a mental model of the field,
  so she accomplishes specific tasks by rote,
  and is unable to compensate for even small differences in tasks.

- A *competent practitioner* is someone whose mental model is robust enough
  to allow her to accomplish normal tasks successfully under normal circumstances.
  Most drivers,
  for example,
  are competent:
  they can handle a car in normal city traffic without worrying about crashing.

- An *expert* is someone who can handle abnormal tasks and abnormal circumstances.
  As we discuss below,
  she doesn't just know more than a competent practitioner:
  her mental model contains many more links between facts.

This model explains the key difference between novice and competent practitioners in any field.
To paraphrase Tolstoy,
everyone who understands something understands it the same way,
but everyone who misunderstands,
misunderstands differently.
The reason is that by definition,
novices don't yet have the right conceptual categories for what they're learning,
so they are putting new knowledge into old boxes,
and those boxes are the wrong ones for this domain.
Everyone who has taught novices has seen this first-hand:
their questions often don't make sense because they've put what they know into the wrong boxes.

As mentioned above,
the key difference between competent practitioners and experts isn't that the latter know more facts:
it's that there are many more linkages *between* facts in an expert's mind.
If we imagine the brain storing knowledge as a graph
(which it emphatically doesn't, but it's a useful metaphor),
then the expert's graph is much more densely connected than that of someone who's "merely competent".
As a result,
the expert is able to jump directly from A to E
when a competent practitioner would have to reason her way through B, C, and D.

This three-stage model has several implications for teaching.
As an example of the first,
our introduction to the Unix shell only introduces 15 commands in two and a half hours.
Ten minutes per command may seem glacially slow,
but that's because our real goal is
to help learners construct the mental model and notional machine that these commands fit into.
That model includes things like:

- If you repeat things manually you'll eventually get them wrong,
  so let the computer repeat things for you by using tab completion and the `history` command.
- Lots of little tools, freely combined,
  are more productive than a handful of "kitchen sink" programs.

As the points above illustrate,
learning consists of more than "just" conveying facts and techniques:
it's just as important to create linkages between concepts.
Telling people that they shouldn't repeat things,
and that they should try to think in terms of little pieces loosely joined,
both set the stage for the introduction of functions in our programming lesson;
explicitly referring back to the shell lesson at that point helps solidify both.

Second,
this model explains why experts are often poor teachers.
Their "knowledge graphs" are so densely connected that they can jump directly to a solution:
it really *is* "obvious" to them.
They then have to struggle to reverse engineer a step-by-step explanation of their thinking
for the novices or competent practitioners they're teaching,
who don't have those links,
and those post hoc rationalizations are often not how people at earlier stages of development
would actually think their way through the problem.

This model also explains why "broadcast mode" teaching is so ineffective.
The transition from novice to competent practitioner
requires the construction of a working mental model of the problem.
Replaying recorded content is as ineffective for most learners as watching TV,
or as a professor standing in front of 400 people in a lecture hall,
because neither can intervene to clear up specific learners' misconceptions.
People who happen to have good enough conceptual categories for a subject
or form them early on
are able to progress through such instruction,
but most others are not.

> #### The Documentation Is Useless
> 
> This is also one of the reasons software documentation is so often frustrating.
> Reference material is opaque to someone who doesn't know what they're looking for,
> such as a novice who doesn't yet have a mental map of the domain.
> On the other hand,
> tutorials meant to help people build such a map
> are too slow and too diffuse for people who already have one.
> It is possible to craft something that serves both communities,
> but it's often simpler to address their needs separately.

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

## Scaling Instruction

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

## Improving Instruction

Teaching is a performing art.
As in music or athletics,
anyone who wants to do it well can draw on an enormous body of experience and research to improve,
but the most important thing is to realize that feedback makes perfect.

Elizabeth Green's book
*[Building a Better Teacher](http://www.amazon.com/Building-Better-Teacher-Teaching-Everyone/dp/0393081591/)*
is a good popular explanation of why.
Its thesis is simple:
most people incorrectly assume that great teachers are born, not made.
From politicians to researchers and teachers themselves,
most of us speak and act as if there's a gene for teaching that someone either has or doesn't.
Most reforms are therefore designed to find and promote those who can and eliminate those who can't.

The problem is,
this assumption is wrong,
so educational reforms based on it (mostly) fail.
Reforms based on changing the culture of teaching would have a greater chance of succeeding,
but as with any cultural change,
they would require the kind of long-term commitment that our society doesn't seem to be very good at.

The book is written as a history of the people who have put that puzzle together in the US,
including Nate Gage and Lee Shulman in the 1960s and 1970s,
Deborah Ball, Magdalene Lampert, and others at Michigan State in the 1980s and 1990s,
and educational entrepreneurs like Doug Lemov today.
Its core begins with a discussion of what James Stigler discovered during a visit to Japan in the early 1990s:

> Some American teachers called their pattern "I, We, You":
> After checking homework,
> teachers announced the day's topic,
> demonstrating a new procedure (I)...
> Then they led the class in trying out a sample problem together (We)...
> Finally, they let students work through similar problems on their own,
> usually by silently making their way through a worksheet (You)...
>
> The Japanese teachers, meanwhile, turned "I, We, You" inside out.
> You might call their version "You, Y'all, We."
> They began not with an > introduction,
> but a single problem that students spent ten or twenty > minutes working through alone (You)...
> While the students worked,
> the teacher wove through the students' desks,
> studying what they came up with and taking notes to remember who had which idea.
> Sometimes the teacher then deployed the students to discuss the problem in small > groups (Y'all).
> Next, the teacher brought them back to the whole group,
> asking students to present their different ideas for how to solve the problem on the chalkboard...
> Finally, the teacher led a discussion,
> guiding students to a shared conclusion (We).

It's tempting to think that this particular teaching technique is Japan's secret sauce:
tempting, but wrong.
The actual key is revealed in the description of Akihiko Takahashi's work.
In 1991,
he visited the United States in a vain attempt to find the classrooms described a decade earlier
in a report by the National Council of Teachers of Mathematics (NCTM).
He couldn't find them.
Instead, he found that American teachers met once a year (if that) to exchange ideas about teaching,
compared to the weekly or even daily meetings he was used to.
What was worse:

> The teachers described lessons they gave and things students said,
> but they did not *see* the practices.
> When it came to observing actual lessons --- watching each other teach ---
> they simply had no opportunity...
> They had, he realized, no *jugyokenkyu*.
> Translated literally as "lesson study",
> *jugyokenkyu* is a bucket of practices that Japanese teachers use to hone their craft,
> from observing each other at work to discussing the lesson afterward
> to studying curriculum materials with colleagues.
> The practice is so pervasive in Japanese schools that it is...effectively invisible.
>
> And here lay the answer to [Akihiko's] puzzle.
> Of course the American teachers' work fell short of the model set by their best thinkers...
> Without *jugyokenkyu*, his own classes would have been equally drab.
> Without *jugyokenkyu*, how could you even teach?

So what does *jugyokenkyu* look like in practice?

> In order to graduate,
> education majors not only had to watch their assigned master teacher work,
> they had to effectively replace him,
> installing themselves in his classroom first as observers and then,
> by the third week, as a wobbly...approximation of the teacher himself.
> It worked like a kind of teaching relay.
> Each trainee took a subject,
> planning five days' worth of lessons... [and then] each took a day.
> To pass the baton,
> you had to teach a day's lesson in every single subject:
> the one you planned and the four you did not...
> and you had to do it right under your master teacher's nose.
> Afterward,
> everyone --- the teacher, the college students, and sometimes even another outside observer ---
> would sit around a formal table to talk about what they saw.
>
> [Trainees] stayed in...class until the students left at 3:00 pm,
> and they didn't leave the school until they'd finished discussing the day's events,
> usually around eight o'clock.
> They talked about what [the master teacher] had done,
> but they spent more time poring over how the students had responded:
> what they wrote in their notes;
> the ideas they came up with, right and wrong;
> the architecture of the group discussion.
> The rest of the night was devoted to planning...
>
> ...By the time he arrived in [the US],
> [Akihiko had] become...famous...
> giving public lessons that attracted hundreds,
> and, in one case, an audience of a thousand.
> He had a seemingly magical effect on children...
> But Akihiko knew he was no virtuoso.
> "It is not only me," he always said...
> "*Many* people."
> After all, it was his mentor...who had taught him the new approach to teaching...
> And [he] had crafted the approach along with the other math teachers
> in [his ward] and beyond.
> Together, the group met regularly to discuss their plans for teaching...
> [At] the end of a discussion,
> they'd invite each other to their classrooms to study the results.
> In retrospect,
> this was the most important lesson:
> not how to give a lesson, but how to study teaching,
> using the cycle of *jugyokenkyu* to put...work under a microscope and improve it.

Putting work under a microscope in order to improve it is commonplace in sports and music.
A professional musician, for example,
will dissect half a dozen different recordings of "Body and Soul" or "Yesterday" before performing it.
They would also expect to get feedback from fellow musicians during practice and after performances.
Many other disciplines work this way too:
the Japanese drew inspiration from [Deming](https://en.wikipedia.org/wiki/W._Edwards_Deming)'s ideas
on continuous improvement in manufacturing,
while the adoption of code review over the last 15 years
has done more to improve everyday programming than any number of books or websites.

But this kind of feedback isn't part of teaching culture in North America.
Here, what happens in the classroom stays in the classroom:
teachers don't watch each other's lessons on a regular basis,
so they can't borrow each other's good ideas.
The result is that *every teacher has to invent teaching on their own*.
They may get lesson plans and assignments from colleagues,
the school board,
a textbook publisher,
or the Internet,
but each teacher has to figure out on their own how to combine that with
the theory they've learned in education school
to deliver an actual lesson in an actual classroom for actual students.

This is a pretty accurate description of what Software Carpentry instructors have to do as well.
Our training course covers the basics of educational psychology and instructional design,
but doesn't walk trainees through our existing material or how to deliver it.
That oversight is completely my fault,
and doubly embarrassing given how many people have asked me to do precisely that.
More importantly,
while instructors may pick ideas up from their fellow teachers at particular bootcamps,
we don't systematically compile and share their experiences.
We don't even really know how much of each topic gets covered in the average bootcamp,
much less how well it works.

Another problem is that there's no longer a "reference implementation" for Software Carpentry.
In January 2013, almost every instructor had taught with Greg Wilson at least one.
That fraction is now closer to a third,
and is inevitably going to continue to drop.
Greg's approach to teaching is almost certainly not the best possible,
but at least it was a *single* way that everyone used to have
as a common starting point.

If we were all at the same university,
or even in the same state or province,
we could try fixing all this with some *jugyokenkyu*.
Unfortunately,
what's described in *Building a Better Teacher* depends on
frequent, time-consuming, high-bandwidth interaction,
but our instructors only teach a couple of times a year as a sideline to their real work.
Outside the days they're actually teaching,
they're scattered across six continents,
trying to catch up on all the work they pushed aside to go run a bootcamp.
