---
layout: lesson
root: ..
title: Collaborative Lesson Development
---
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
