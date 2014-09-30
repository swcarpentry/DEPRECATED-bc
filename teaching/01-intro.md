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

> ### What We Leave Out
> 
> This training takes a narrow "cognitivist" perspective on education.
> We can afford to do that because most of our learners are privileged:
> well-fed, physically safe, already well-educated, highly motivated, etc.
> There are many other perspectives ---
> see [this site](http://www.learning-theories.com/) for summaries ---
> and it's worth exploring them for insights like this:
>
> "If poor inner-city children consistently outscored children
> from wealthy suburban homes on standardized tests,
> is anyone naive enough to believe that we would still insist on using these tests
>  as indicators of success?"
> 
> --- Kenneth Wesson, in Littky and Grabelle's *The Big Picture*

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

Here's an example.
Today,
most children are taught to read using a bottom-up technique called
[phonics](http://en.wikipedia.org/wiki/Phonics),
which introduces them to the sounds of letters,
then puts those letters and their sounds together to make simple words like "bug" and "run",
then assembles the words into sentences and stories.
But another technique called [whole language](http://en.wikipedia.org/wiki/Whole_language)
takes a different approach.
Instead of starting with letters,
its proponents teach children to recognize entire words
so that they get the positive feedback of being able to read whole stories as early as possible.
The basic idea is that once children find reading rewarding,
they'll be motivated to go back and figure out the letters.

Both approaches are consistent with what we know about how children learn.
To determine which is better,
we need to study them in action.
The findings have been contradictory:
where one series of studies will find that phonics works best,
another will find that whole language yields better results.

One explanation for these confusing results is that a third factor is at play.
If children sense that their teachers are enthusiastic about something,
they will respond to that enthusiasm and learn no matter what technique is being employed.
If they feel their teachers are just marking time,
though,
they will emulate that instead.
Thus,
when studies are done by people who are proponents of a new technique
(as they often are),
they are likely to produce positive results,
but when those techniques are deployed more widely
to teachers who have seen one bandwagon after another roll by,
results will regress to the mean.

> #### The Song Remains the Same
>
> This same debate between advocates of a bottom-up "principles first" approach
> and people who emphasize early rewards and motivation
> plays out in programming education as well.
> Software Carpentry takes the second approach
> because we can't start with foundational theory
> and get novices anywhere useful in a two-day workshop.
> Our position is bolstered by research by people like Mark Guzdial and Barbara Ericson
> on media-first computation,
> which we will discuss later.

## The Most Important Number in Teaching (and Programming)

The next piece of psychological theory we need concerns memory.
As a very rough approximation,
human memory is divided into two layers.
The larger is called *long-term memory* or *persistent memory*,
and its most important characteristics are:

1.  It is essentially *unbounded*: a healthy adult will die before her long-term memory fills up.
2.  It is *associative*, i.e., it works by pattern-matching.
3.  It is *slow* --- too slow to support conscious thought and action.

To accommodate the third of these facts,
evolution has put a second layer on top of long-term memory
called *short-term* or *working* memory.
It is also associative,
but it's much faster.
It is also much, much smaller:
[early estimates](https://en.wikipedia.org/wiki/The_Magical_Number_Seven,_Plus_or_Minus_Two) were that
it could hold 7&plusmn;2 items at once,
but more recent research puts the number even lower.

One sign of how small working memory is can be found in the length of phone numbers.
Back when phones had rotary dials,
people had to be able to read a string of more-or-less random digits in a phone book
and hold them in their head long enough for the phone's dial to go around several times.
At six, seven, or eight digits,
most adults could complete the task without any trouble;
at twelve digits,
most would forget the last couple of digits by the time they needed to dial them in.

If this is true,
though,
how are we still able to dial numbers
when most of them include an area code?
The answer lies in another phenomenon called [chunking](https://en.wikipedia.org/wiki/Chunking_%28psychology%29).
If we see several things together often enough,
our brains will group them into a chunk,
then remember the chunk as one thing.
For example,
most of us don't remember our names as a sequence of letters like 'S', 'e', 'v', 'e', 'r', 'u', 's'.
Instead,
we have seen it so often that it becomes one thing: 'Severus'.
This allows us to work with much more information --- groups of groups of groups --- but
there is always the risk that our brains will "see" chunks that aren't actually there.
For example,
if you flash a playing card with this pattern on it:

~~~
*   *
  *
  *
*   *
~~~

most people will remember seeing five spots,
rather than six,
because their brain has learned that the 'X' pattern contains five pieces.

> #### Helping People Form Chunks
>
> Telling people, "This is a pattern," isn't enough to form a chunk:
> they need to see it repeatedly for the pattern to sink in.
> What we *can* do is show them lots of examples
> that contain patterns we want them to learn to recognize,
> like:
>
> ~~~
> for item in collection:
>     if item passes test:
>         temp = do something to item
>         result.append(temp)
> ~~~
>
> In order for the chunk to form,
> we need to present it in the same way every time:
> an experienced programmer will recognize that the code above is equivalent to:
>
> ~~~
> for item in collection:
>     if item passes test:
>         result.append(do something to item)
> ~~~
>
> but a novice's pattern-forming will be disrupted by the "superficial" difference.

Another sign of our brains' limitations is the fact that
all sports teams have about half a dozen members.
Larger "teams" alwyas break down into sub-groups:
the forwards and backs in rugby,
the infield and outfield in baseball,
and so on.
A squad with half a dozen members has also been
the basic unit of pretty much every military organization
throughout recorded history
because our short-term memory can't keep track of any more people than that.

The main implication of this is that we should teach in very small chunks.
Ideally,
learners should only have to digest a handful of new ideas at once.
Those ideas should immediately be reinforced through practice
to help transfer them from short-term to long-term memory.
This is very different from the hour-long lecture format
that is commonly used for "intellectual" subjects,
but is instantly familiar from music lessons and other "performance" subjects.
The former may give instructors (the illusion of) more control,
but the latter produces better results.

> #### Live Coding
>
> Software Carpentry workshops don't use slides to teach programming.
> Instead,
> we use live coding,
> i.e.,
> we share our screen and type code in right in front of the class.
> This has several benefits:
>
> 1. *It encourages us to break things up into smaller chunks.*
>    If learners are typing along with instructors,
>    it's easy for the latter to say "now try this" at frequent intervals.
>
> 2. *It slows us down.*
>    It's all too easy to race ahead of learners when flipping through slides;
>    if the instructor has to type things in herself,
>    she's less likely to do that.
>
> 3. *Learners get to see how we use programming tools.*
>    It's hard for a slide to show tab completion
>    or the way we locate a function's definition.
>
> 4. *Learners get to see how we recover from mistakes.*
>    We're often told that one of the most useful parts of a workshop is
>    seeing how instructors get out of the holes they've dug for themselves,
>    since static tutorials (in books, slides, or online)
>    usually show only the "happy path".
>
> 5. *It gives learners permission to make mistakes themselves.*
>    If learners see an instructor making mistakes,
>    they're likely to be more forviging of their own,
>    and less likely to become discouraged.

## Stages of Cognitive Development

The next piece of background we need is
how learners' mastery of a subject changes.
In the 1980s,
Patricia Brenner developed
[a five-stage model of skill acquisition](http://www.amazon.com/From-Novice-Expert-Excellence-Commemorative/dp/0130325228/)
that has been applied successfully in many fields.
We will simplify that model to three stages:

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

## Implications

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
This leads to "expert blind spot":
experts have forgotten what it's like to not understand something,
so they can't "see" the problems that novices are having.

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
