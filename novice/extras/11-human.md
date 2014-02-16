---
layout: lesson
root: ../..
title: The Human Side of Things
level: novice
---
In my experience, the things that go wrong most often in software
development projects&mdash;undergraduate or professional&mdash;have nothing to
do with software.  Instead, the worst mistakes people make are all
related to the human side of things.  Since some of you will only read
the first few pages of any document, I have put the three most
important topics here in a chapter of their own.

#### Meetings

I have a theory.  It's not about why there are no toilets on the
*Enterprise* (hint: transporters).  It's about why bad
software exists.

> Fact #1: the two biggest reasons software projects blow up are poor
> requirements management and integration failure: either the wrong
> thing is built, or the pieces don't work together when they're
> assembled.
>  
> Fact #2: the usual way to sort out requirements and negotiate
> interfaces is face-to-face discussion.  That means having meetings.
> 
> Fact #3: most people are really bad at meetings: they don't have an
> agenda going in, they don't take minutes, they waffle on or wander
> off into irrelevancies, they repeat what others have said or recite
> banalities simply so that they'll have said something, they hold
> side conversations (which pretty much guarantees that the meeting
> will be a waste of time), and so on.
> 
> Conclusion: people create bad software because they're bad at
> meetings.
> 
> Corollary: improve people's meeting skills, and you'll improve the
> software they produce.

So, here's how to run a meeting.

*Start with an agenda.* If nobody cares enough about the meeting
to write a point-form list of what's supposed to be discussed, the
meeting itself probably doesn't need to happen.  Agendas also help you
keep the meeting moving (as in, "That's very interesting, Yousuf, but
we have five other topics to get through in the next fifteen minutes,
so could you please make your point?").

*Only one person talks at a time.*  That means no side
conversations, no catching up with blogs, and no cell phones or
Blackberries.  This rule is important because:

*   It's more efficient.  Nobody can actually pay attention to two
    things at once: if distractions are tolerated, people will miss
    things, or they'll have to be repeated, or both.

*   Not paying attention is insulting&mdash;chatting with a friend on
    IM while someone explains what they did last week is a really
    effective way to say, "I don't think your work is important."

*Put someone in charge of the meeting.* "In charge" means
keeping the meeting moving, glaring at people who are muttering to one
another or checking email, and telling people who are talking too much
to get to the point.  It does *not* mean doing all the talking;
in fact, whoever is in charge will usually talk less than anyone else,
just as a referee usually kicks the ball less often than the players.

*No one is allowed to interrupt* (except the person in charge).
Raise a finger (no, not that one), or make one of the other gestures
people make at high-priced auctions instead if you want to speak next.
If the speaker doesn't notice you, the person in charge ought to.  Oh,
and while other people are talking, take notes of questions you want
to ask or points you want to make.  You'll be surprised how smart it
makes you look when it *is* your turn to speak.

*Have one person take minutes,* and circulate them (by email, or
in a blog or a wiki) by the end of the day of the meeting.  These
minutes aren't supposed to be a word-for-word record of who said what.
Instead, they should summarize the main points people made, the
questions that no one could answer, and the things people promised to
do, so that:

*   *People who weren't at the meeting can keep track of
    what's going on.*  You and your fellow students all have to juggle
    assignments from several other courses while doing this project,
    which means that sometimes you won't be able to make it to team
    meetings.  A wiki page, email message, or blog entry is a much more
    efficient way to catch up after a missed meeting or two than asking
    a team mate, "Hey, what did I miss?"

*   *Everyone can check what was actually said or promised.*
    More than once, I've looked over the minutes of a meeting I was in
    and thought, "Did I say that?" or, "Wait a minute, I didn't
    promise to have it ready then!"  Accidentally or not, people will
    often remember things differently; writing it down gives team
    members a chance to correct mistaken or malicious interpretations,
    which can save a lot of anguish later on.

*   *People can be held accountable at subsequent meetings.*
    There's no point making lists of questions and
    action items if you don't follow up on them later.  If you're
    using a ticketing system, the best thing to do is to create a
    ticket for each new question or task right after the meeting,
    and update those that are being carried forward.  That way,
    your agenda for the next meeting can start by rattling through
    a list of tickets.

That's it.  Run all your meetings like this for a month, with the goal
of making each one a minute shorter than the one before, and I promise
that you'll build better software.

Oh, wait, there's one more rule: don't meet just to exchange
information.  Once you're using version control and tickets
you'll be able to keep track of your teammates' actual progress
via the web.  Meetings should be held to make decisions about
design, work allocation, and what to cut when time gets tight.
Remember, you can read faster than anyone can speak: if someone
has facts for the rest of the team to absorb, the most polite
way to communicate them is to type them in.

#### Crunch Mode

I used to brag about the hours I was working.  Not in so many words,
of course&mdash;I had *some* social skills.  Instead, I'd show up for
class around noon, unshaven and yawning, and casually mention how I'd
been up 'til 6:00 a.m. hacking away at some monster bug or other.

Looking back, I can't remember who I was trying to impress.  Instead,
what I remember is how much of the code I wrote in those all-nighters
I threw away once I'd had some sleep, and how much damage the bugs I
created in those bleary-eyed stupors did to my grades.

My mistake was to confuse "working" with "being productive".  You
can't produce software (or anything else) without doing some work, but
you can easily do lots of work without producing anything of value.
Scientific study of the issue goes back to at least the 1890s.  The
most important results for students are:

*   Working more than eight hours a day for an extended period of
    time lowers your total productivity, not just your hourly
    productivity&mdash;i.e., you get less done in total (not just per hour)
    when you're in crunch mode than you do when you work regular hours.

*   Working over 21 hours in a stretch increases the odds of you
    making a catastrophic error just as much as being legally drunk.

These facts have been reproduced and verified through hundreds of
experiments over the course of more than a century.  The data behind
them is as solid as the data linking smoking to lung cancer.  However,
while most smokers will admit that their habit is killing them, people
in the software industry still talk and act as if they were somehow
exempt from these findings.  To quote Robinson's article:

> When Henry Ford famously adopted a 40-hour workweek in 1926, he was
> bitterly criticized by members of the National Association of
> Manufacturers. But his experiments, which he'd been conducting for
> at least 12 years, showed him clearly that cutting the workday from
> ten hours to eight hours&mdash;and the workweek from six days to five
> days&mdash;increased total worker output and reduced production
> cost... the core of his argument was that reduced shift length
> meant more output.
>  
> ...many studies, conducted by businesses, universities,
> industry associations and the military, ...support the basic
> notion that, for most people, eight hours a day, five days per week,
> is the best sustainable long-term balance point between output and
> exhaustion. Throughout the 30s, 40s, and 50s, these studies were
> apparently conducted by the hundreds; and by the 1960s, the benefits
> of the 40-hour week were accepted almost beyond question in
> corporate America. In 1962, the Chamber of Commerce even published a
> pamphlet extolling the productivity gains of reduced hours.
>  
> But, somehow, Silicon Valley didn't get the memo...

I was part of a data visualization startup in the mid-1990s.  Three
months before our first release, the head of development "asked" us to
start coming in on Saturdays.  We were already pulling one late night
a week at that point, and most of us were also working at least a
couple of hours at home in the evenings.  It's hardly surprising that
we missed our "can't miss" deadline by ten weeks, and had to follow up
our 1.0 release with a 1.1, and then a 1.2, in order to patch the
crash-and-lose-data bugs we'd created.  We were all zombies, and
zombies can't code.

Those kinds of hours are sadly still normal in many parts of the
software industry, and also in university programs.  Everyone knows
that designing and building software is a creative act that requires a
clear head, but many of those same people then act as if it was like
digging a ditch.

The big difference is that it's hard to lose ground when digging
(though not impossible).  In software, on the other hand, it's very
easy to go backward.  It only takes me a couple of minutes to
create a bug that will take hours to track down later&mdash;or days, if
someone else is unlucky enough to have to track it down.  This is
summarized in Robinson's first rule:

> Productivity varies over the course of the workday, with the
> greatest productivity occurring in the first four to six
> hours. After enough hours, productivity approaches zero; eventually
> it becomes negative.

It's hard to quantify the productivity of programmers, testers, and
UI designers, but five eight-hour days per week has been proven to
maximize long-term total output in every industry that has ever been
studied.  There's no reason to believe that software development is
any different, or that student programming is different from full-time
programming in industry.

Ah, you say, that's "long-term total output".  What about short
bursts now and then, like pulling an all-nighter to meet a deadline?
Well, that's been studied too, and the results aren't pleasant.  Your
ability to think drops by 25% for each 24 hours you're awake.  Put it
another way, the average person's IQ is only 75 after one all-nighter,
which puts them in the bottom 5% of the population.  Two all nighters
in a row, and their effective IQ is 50, the level at which people are
usually judged incapable of independent living.

The catch in all of this is that
*people usually don't notice their abilities declining*.
Just like drunks who think they're still able to drive, people who are
deprived of sleep don't realize that they're not finishing their
sentences (or thoughts).  They certainly don't realize that they're
passing parameters into function calls the wrong way around, or that
what they're typing in will all have to be deleted and re-done
tomorrow, when it will take longer than it would have if they'd just
gone home and gotten a good night's sleep.

The moral of this story?  Think very hard about what's more important
to you: the amount of good work you produce, or how much of a martyr
you appear to be.  Then think about which of those other people are
actually going to care about, and pace yourself accordingly.

#### Time Management

"But&mdash;but&mdash;I have so many assignments to do!", you say.  "And
they're all due at once!  I *have* to work extra hours to get
them all done!"

No.  You can work *smarter* instead.  I'm less intelligent than
many of the people I've worked with over the years.  I still manage to
get a lot done because I learned early on that organization and focus
are more important than raw IQ.

In order to be productive, you have to do two things: prioritize, and
focus.  The first is important because people are naturally very good
at spending hours on things that don't need to be done, and then
finding themselves with too little time for the things that actually
count.  It can actually be expressed as an algorithm:

*   Make a list of the things you have to do.  I still use a
    hardcover lab notebook for this, since I can doodle in it when
    I'm on the subway, but a lot of people keep a personal wiki,
    or send themselves email messages that then go into a folder
    titled "To Do".  However you do it, the important thing is
    to *write it all down*.  Your mind can only keep seven
    or so things in your short-term memory at once; if you try to
    manage a to-do list longer than that in your head, you will
    forget things.

*   Weed out everything that you don't need to do right away.
    Notice that I said "need", not "want": if you want to mess
    around with a new blogging tool, that's fine, but that's play time,
    not work time, and we're talking about getting work done.

*   Prioritize what's left by sorting the list so that the most
    important tasks are at the top.  (I don't worry about getting the
    stuff below the first three or four lines into exact order, since
    I'm going to re-check my list before I get to them anyway.)

*   Make sure you have everything you need to see the first task
    through: files from version control, the assignment
    specification, whatever libraries and documentation you need,
    a fresh cup of coffee, a comfortable chair, etc.  Don't give
    yourself an excuse to interrupt your own work: the world will
    provide enough of those.

*   Shut down your mail client, and turn off instant messaging and
    your cell phone.  Don't panic, it's only for an hour&mdash;most people
    can't stay focused longer than that, and anyway, you'll need to
    stretch your muscles and get rid of that coffee you drank.

*   Set an alarm to go off in sixty minutes, and *focus*.
    Don't switch tasks in that hour unless you absolutely have to.
    Instead, if you remember an email message that needs to be sent, or
    discover a couple of tests that really should be written, add a note
    to your to-do list.  (This is another reason I keep mine in a lab
    notebook: the few seconds it takes to pick up a pen and jot
    something down gives my hands a rest from the keyboard.)

*   When your hour is up, take a break: check mail (but don't
    reply to anything that isn't urgent), go to the washroom, stretch a
    little, and then re-order your to-do list and start the next round.

If any task on your list is more than an hour long, break it down into
into smaller pieces and prioritize those separately.  Keep in mind
that the future is approaching at a fixed rate of one day every 24
hours: if something's going to take sixty hours to do, you'd better
allow at least ten working days for it, which means you'd better tackle the
first piece two working weeks before the deadline.  And since breaking
large tasks down into small ones takes time, don't be surprised if
"plan XYZ" appears as a task in your list.

The point of all this organization and preparation is to get
yourself into the most productive mental state possible.
Psychologists call it flow; athletes call it "being in the
zone", while musicians talk about losing themselves in what
they're playing.  Whatever name you use, you will produce much
more per unit of time in this state than normal.

That's the good news.  The bad news is that it takes roughly ten
minutes to get back into a state of flow after an interruption, no
matter how short the interruption was.  This means that if you are
interrupted half a dozen times per hour, you are *never* at
your productive peak.  It's very much like processes being paged in
and out by an operating system: if it happens too often, the CPU
spends all its time moving things around, and none doing useful
work.

As a concrete example, a student of mine kept a stopwatch beside his
computer for a couple of weeks during term.  Every time he read mail,
put his instant messaging client in the foreground, or went to
Manchester United's web site, he hit the button to stop it.  At the
end of two weeks, he discovered that he only spent 28% of his
"working" time working.  Put it another way, he could have finished
his assignments in a third of the time they actually took.

If you want to get as much work done as possible, so that you will
have enough time for other things, it is therefore crucial that you
separate the two.  IM and your cell phone are great for socializing;
it's easy to tell yourself that you should leave them on so that your
teammates can reach you, but it just ain't true.  The best way to help
your teammates is to get your part of the project done; the best way
to show your significant other that you really love them is to finish
your assignments as quickly as you can so that you can go out and see
the latest Pixar movie without having to cut out and go back to the
lab.

Making lists and setting one-hour alarms will probably seem a little
earnest at first.  Trust me, though: your friends will stop mocking
you once they see that you're able to finish your assignments and
still have time to play some badminton and catch a movie.  They may
even start to imitate you.

#### Are You Sitting Comfortably?

A complete working environment needs more than just software.
Unfortunately, most university labs seemed designed to make everything
below difficult or impossible to achieve.

*   *Peace and quiet.* Study after study has proved that this has more
    impact on productivity than a fast network, a fat disk, or
    caffeine, but most workplaces are still too crowded, too noisy,
    and filled with too many interrupts.  As mentioned above, it takes
    most people ten minutes to get back into a state of flow after
    being distracted, which means that half a dozen interruptions per
    hour effectively renders someone zero percent effective.  I know
    people say, "If I can't overhear what other people are talking
    about, I might miss something important," but that only applies if
    the only people you're overhearing are members of your own team
    (and even then, it's a dubious claim).

*   *Comfortable seating.* A good chair with a firm back
    costs one fifth of what a mid-range PC does.  A full-sized keyboard
    (I have large hands&mdash;most laptop keyboards force me to bend my
    wrists uncomfortably) costs fifty dollars, and a lamp with a
    soft-light bulb is another forty.  The combination doesn't just let
    me program longer each day; it also helps ensure that I'll still be
    able to program five or ten years from now without chronic back and
    wrist pain.  Compare this with what's in most computer labs: their
    lighting gives glare without illumination, the dark desktops make
    the optical mice jerky, and the low-budget chairs are guaranteed to
    make your lower back ache after an hour.

*   *A pad of gridded paper and several ballpoint pens.*  I
    often make notes for myself when programming, or draw box-and-arrow
    diagrams of my data structures when debugging.  I used to keep an
    editor open in a background window to do the former, but when my
    wrists started acting up, I discovered that taking my hands away
    from the keyboard for a few moments to scribble something down
    provided welcome relief.  I also quickly discovered that the odds of
    me being able to read my own notes the next day rose dramatically if
    I used gridded paper to line them up.

*   *A heavy mug for coffee or tea.* I don't know why, but a
    styrofoam cup, or a normal teacup, just isn't as satisfying as a
    little hand-sized ceramic boulder.  Maybe it satisfies my
    subconscious Neanderthal urge to club my computer to death when it
    misbehaves...

*   *A rubber duck.* One of computing's many apocryphal
    stories holds that someone&mdash;Brian Kernighan, maybe, or Dennis
    Ritchie&mdash;keeps a rubber duck next to his computer.  Whenever a bug
    takes more than a few minutes to track down, he puts the duck on his
    desk and explains the problem to it.  Why?  Because speaking out
    loud forces you to marshal your thoughts, which in turn highlights
    any contradictions or missed steps that you hadn't noticed while
    everything was just swirling around inside your head.

*   *A squirt bottle of glass cleaner and a box of kleenex.*
    I can't stand smears on my screen.  They drive me nuts.  Whenever
    I'm showing something to someone, and they actually *touch my
    screen* instead of just pointing, I have another Neanderthal
    fantasy, except this time it's not subconscious...

*   *A chess set*, or a deck of cards, or some juggling
    balls.  I'm a very bad chess player.  Luckily, so are most
    people, so it's usually possible to find someone at my level
    for a ten-minute game at lunch.  Other programmers I know play
    euchre, or knit&mdash;a programmers' stitch 'n' bitch session
    can be jaw-dropping to listen to.  Few people can focus for
    more than a few hours before their productivity drops; it's
    better to acknowledge this, and take a break in the middle of
    the day, than to say, "Must... keep... coding..." and produce
    garbage that just has to be rewritten later.

*   *Books.* You can find a lot on-line, but it's hard to
    google with your feet up on your desk, and even harder to fold down
    the corners on web pages.  When I want API documentation, I use the
    web; when I want a tutorial, I still prefer the printed page.

*   *Running shoes.* Back when I was a part-time grad
    student, I had a settled routine: I brought three sets of gym gear
    to the office at the start of the week, worked out at lunchtime on
    Monday, Wednesday, and Friday, and took my stuff home at the end of
    the week.  After two months of this, I came in to find that my
    co-workers had hung a little chandelier made of air fresheners over
    my desk.  Since then, I've rented a locker at the gym, but I still
    try to get some exercise several times a week&mdash;it helps my
    concentration and stamina a lot more than any amount of coffee..

*   *Pictures.* Ah, the nesting instinct.  Everyone wants to
    feel at home; everyone wants to make wherever they are uniquely
    theirs.  I hang a few postcards on the wall wherever I work, along
    with a photograph of my wife and daughter taken ten hours after she
    was born (my daughter, that is), just to remind me what's really
    important.
