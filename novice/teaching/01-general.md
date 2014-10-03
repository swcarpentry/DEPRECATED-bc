---
layout: lesson
root: ../..
title: General Advice
level: novice
---
This is a placeholder for general notes about teaching.
For up-to-date information about software installation and configuration problems,
and their solutions,
see [this wiki page](https://github.com/swcarpentry/bc/wiki/Configuration-Problems-and-Solutions).

#### Teaching Notes

*   For workshops that extend over more than two days
    (e.g., four afternoons spread over two weeks),
    it's a good idea to email the learners at the end of each day with a summary of what was taught
    (with links to the relevant online notes).
    This allows absent learners to catch up before the next session,
    and is also a great opporunity to present the lessons of the day
    in the context of the entire workshop.

*   Point learners at [http://software-carpentry.org/v5/](http://software-carpentry.org/v5/),
    which is the permanent home of the current learning materials,
    and at [http://software-carpentry.org/v4/](http://software-carpentry.org/v4/),
    which is where our previous materials live.
    The former corresponds to what they're being taught;
    the latter covers more ground in video as well as in slides and prose.
    They should also be direct to
    [Software Carpentry's FAQ](http://software-carpentry.org/faq.html).

*   Explain that the lesson materials can all be freely re-mixed and re-used
    under the [Creative Commons - Attribution](../../LICENSE.html) (CC-BY) license,
    provided people cite us as the original source
    (e.g., provide a link back to our site).
    However,
    Software Carpentry's name and logo are trademarked,
    and cannot be used without our permission.
    We normally grant this to any class that
    (a) covers our core topics and
    (b) has at least one badged instructor on the teaching roster,
    but are happy to discuss specifics.

*   Plan for the first 30-60 minutes of the workshop to be spent on installation and setup,
    because it's going to happen anyway.
    Running a pre-workshop "help desk" doesn't really affect this:
    the people who are most likely to have installation problems
    probably won't show up.
    (We fantasize occasionally about turning people away if they haven't installed software,
    or at least downloaded the installers,
    but in practice it's hard to do.)

*   Emphasise that good software development skills contribute to productive, reproducible, reusable research.

*   Have learners post a red sticky note
    on their laptop
    whenever they have a question or need help.
    Have them take down their sticky notes at the start of each practical exercise,
    and then post a green one when they're done
    (and a red one when they need help).

*   At lunch and again at the end of the day,
    ask learners to write one good point (i.e., something they learned or enjoyed)
    on their green sticky note
    and one bad point (i.e., something that didn't work, that they didn't understand, or that they already knew)
    on their red one.
    It only takes a couple of minutes to sort through these,
    and it's a quick way to find out how things are actually going.

*   At the very end of the workshop,
    ask learners to alternately give one good point or one bad one
    about the entire workshop. 
    Write these down on the whiteboard as they come in,
    and do not allow repeats
    (i.e., every point has to be new).
    The first few negative points will usually be pretty innocuous;
    after those have been exhausted,
    you will start to get the real feedback.

*   As a variation on the red/green sticky notes,
    make little name tents out of red and green paper,
    held together with name tag labels.
    The learners write their names on the name tags,
    and prop the tents either green side up or red side up
    depending on the feedback they want to give about the lesson being too fast or too slow.

*   Back up the material with your own anecdotes, experiences and evidence&mdash;it
    makes you more credible,
    helps learners understand how to apply what you're teaching to their own problems,
    and prevents the lectures from becoming too dry.

*   Keep a running list of the commands encountered so far in the lesson
    in the Etherpad
    or on a whiteboard adjacent to the projection screen.
    Encourage learners (particularly ones who already know the material and might otherwise get bored)
    to take notes in the Etherpad as well.
    This reduces the effort per learner,
    gives you a chance to see what they think you're saying,
    and provides a record after the workshop of what was actually taught.

*   When the co-instructor isn't teaching,
    she can answer questions on the Etherpad
    and update it with the key points made by the instructor
    (along with commands
    and any related points the instructor may not have mentioned).
    It's less disruptive to the "live" instructor than interjecting with these points,
    but allows the attendees to get the shared expertise from both instructors.

*   For workshops that extend over more than two days (e.g. four afternoons spread over two weeks),
    it's a good idea to email the learners at the end of each day with a summary of what was taught
    (with links to the relevant online notes).
    Not only does this allow absent learners to catch up before the next session,
    it's also a great opporunity to present the lessons of the day in the context of the entire workshop.

*   The long-form notes are intended as a script for instructors
    and as self-study material for learners.
    Do *not* show these notebooks to learners:
    instead,
    start with a blank notebook when teaching and add code as you go.
    This helps prevent you from racing ahead of learners
    and makes it easier for you to improvise in response to their questions.

*   Point learners at the online versions of the long-form notes
    (either on your workshop's home page
    or at [http://software-carpentry.org/v5/](http://software-carpentry.org/v5/)
    *after* the lesson is done:
    if you do it before the lesson,
    they'll try to read the notes while you're trying to talk.

*   If you're really keen,
    keep the SVG's of the diagrams handy in the directory where you're doing your teaching
    so that you can include them in your notebooks by adding an `<img src="...">` element to a Markdown cell
    (or just display them in your browser).
    Most people don't ever actually do this though,
    either because they forget to
    or because they have a whiteboard or flipchart handy.

*   There are (at least) three ways to get data files to learners at the start of a lesson:
    1.  Create a zip file, add it to your workshop's repository, and put a link to it in your workshop's `index.html` page
        so that they can click, download, and unzip.
        This uses something everyone already understands,
        but does assume they know how to navigate from their download directory to their working (lesson) directory,
        which is often not the case.
    2.  Create a throwaway Git repository on GitHub and tell them to run one command to clone it at the start of class.
        This (usually) works even if they've never used Git,
        and as a bonus,
        lets you identify people who (are going to) have Git problems early.
    3.  Paste the data into an Etherpad for learners to copy.
        As a bonus,
        this lets you identify people who (are going to) have trouble using a text editor early.

*   If you are using multiple windows 
    (e.g., a command window and an editor window)
    make sure they are both large enough to be visible by all attendees.
    Remember to pause when switching from one window to the other
    so that learners don't become confused.
    If possible,
    use different background colors for different text windows
    to make it easier for learners to tell them apart
    (but keep in mind red-green and blue-yellow color blindness).

*   As you type at the command line,
    read out what you're typing.
    Remember that most learners can only go half as fast as you,
    because they have to watch you type
    then type it all in again themselves.

#### Pitfalls

Instructing at a workshop isn't trivial.
The most important thing is to remember that
no lesson plan survives contact with the audience.
Whether it's the network going down
or the sudden realization that many of your learners *don't* know how to use SSH,
you will frequently need to improvise.
And even when there aren't hiccups like those,
try your best to adjust your pace and examples based on
learners' questions, puzzled looks, and sighs of impatience.

**Allow enough time for setup.**

In almost all cases,
learners want to use their own laptops during workshops
so that they leave with a working environment set up.
Even if you ask attendees to prepare beforehand,
and give them detailed instructions,
some will not have time, 
or will have run into problems that they're not yet able to fix.
You should therefore schedule at least 20 minutes for
<em>checking the learners' machines</em>
at the beginning of the first day.
Some workshops start early on the first day to allow time to
troubleshoot setup problems.

**Don't ignore your learners.**

You're not there to reproduce one of our online videos in person:
you're there to interact with people so that they get a
better learning experience.
You shouldn't ever go more than two or three minutes without
asking a question (and listening to the answer),
and if it has been 15 minutes since any of your learners asked one,
odds are you've either lost them or are boring them.

**Don't bore your learners.**

Your audience will never care more about what you're teaching 
than you appear to,
so if they get the feeling you're not interested in it,
they won't be either.
This does <em>not</em> mean you have to shout,
crack three jokes a minute,
or harangue them about how this stuff is really, really important,
but you do owe it to your audience to show up mentally as well as physically.
  
**Don't be all talk, no action.**  

The more time folks spend with their hands on the keyboards 
doing exercises,
the more time they're actually paying attention.
The students have their computers in front of them:
if you talk for more than five minutes without asking them to 
use their computers, they'll do so anyway&mdash;on Facebook.
  
**Don't use magic.**  

Typing too fast, using shortcuts or commands 
learners haven't seen yet&mdash;basically,
any time you say, "Don't worry about this just now,"
or they say, "Wait, how did you do that?" or, 
"Can you please slow down, I can't keep up,"
you're no longer actually helping them.

**Don't take over the keyboard.**

When you go around helping people out as they work
on the material, never put your hands on a learner's
computer.  It's too easy to step in and do a few
quick changes that they can't follow and certain
groups of learners will never say that this makes them
uncomfortable - or ask questions.

**Never copy and paste.**

It's easy to find code snippets on the web
but impress upon learners that part of what
they are learning is the muscle memory of 
doing the text input.  To drive this home, you
could show them an example of `rm -rf /` in a sandbox
just to make it clear that there are risks
to running code you don't understand.

**Don't ignore feedback.**

The feedback you get from learners on sticky notes or through
surveys is pointless if you don't pay attention to it 
(or worse, if you explain it away).
There's no point collecting feedback 
during and after each workshop 
if you don't change what and how you teach to reflect it.

**Tell learners "why".**  

Most of our learners are graduate students in science and engineering,
so they know what evidence looks like, 
and why working practices should be evidence-based.
That doesn't mean you have to have the whole of 
empirical software engineering at your fingertips,
but please do read
*[Facts and Fallacies of Software Engineering](http://www.amazon.com/Facts-Fallacies-Software-Engineering-Robert/dp/0321117425/)*
and sprinkle a few of the findings it quotes into your lessons.
  
**Don't show them the forest but not the trees.**  

The things we teach reinforce each other,
so tie them together at every opportunity.
Point out that connecting things with pipes in the shell is like chaining functions together,
or that they can use a shell script to re-run a bunch of different tests before committing to version control,
and so on.
If possible, take 15 minutes or so each day to show them how you use these tools in your day-to-day work.

**Don't underestimate setup requirements.**  

Do you have enough power outlets?
(Are you sure?)
Do you have enough bandwidth to handle fifty people hitting your version control repository at the same time?
(How do you know?)
Can everyone actually log in?
Are the washrooms unlocked?
Does campus security know you're using the room over the weekend?

**Don't let your learners ignore each other.**

Software Carpentry workshops are a great networking opportunity
for our learners (and for us, too).
Get to know your learners by name,
have them work in pairs, and get them to mix up the pairs
at least a couple of times.
Encourage them to chat to one another at coffee breaks and lunch,
and to get a pizza or some curry together for dinner on the first day.

**Relax.**

Something always fails to install for someone
(or they fail to install anything at all),
or a bunch of learners are accidentally locked out of the building after lunch,
or whoever was supposed to drop off power bars didn't.
Roll with it,
and remember to laugh
(even if it's a bit hysterically).
