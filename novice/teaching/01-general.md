---
layout: lesson
root: ../..
title: General Advice
level: novice
---
This is a placeholder for general notes about teaching.

#### Teaching Notes

*   For bootcamps that extend over more than two days (e.g. four afternoons spread over two weeks), it's a good idea to email the learners at the end of each day with a summary of what was taught (with links to the relevant online notes). Not only does this allow absent learners to catch up before the next session, it's also a great opporunity to present the lessons of the day in the context of the entire bootcamp.

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

*   Plan for the first 30-60 minutes of the bootcamp to be spent on installation and setup,
    because it's going to happen anyway.
    Running a pre-bootcamp "help desk" doesn't really affect this:
    the people who are most likely to have installation problems
    probably won't show up.
    (We fantasize occasionally about turning people away if they haven't installed software,
    or at least downloaded the installers,
    but in practice it's hard to do.)

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

*   At the very end of the bootcamp,
    ask learners to alternately give one good point or one bad one
    about the entire bootcamp. 
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

*   Keep a running list of the commands encountered so far in the lesson
    in the Etherpad
    or on a whiteboard adjacent to the projection screen.

*   When the co-instructor isn't teaching,
    she can answer questions on the Etherpad
    and update it with the key points made by the instructor
    (along with commands
    and any related points the instructor may not have mentioned).
    It's less disruptive to the "live" instructor than interjecting with these points,
    but allows the attendees to get the shared expertise from both instructors.

*   For bootcamps that extend over more than two days (e.g. four afternoons spread over two weeks),
    it's a good idea to email the learners at the end of each day with a summary of what was taught
    (with links to the relevant online notes).
    Not only does this allow absent learners to catch up before the next session,
    it's also a great opporunity to present the lessons of the day in the context of the entire bootcamp.

*   The long-form notes are intended as a script for instructors
    and as self-study material for learners.
    Do *not* show these notebooks to learners:
    instead,
    start with a blank notebook when teaching and add code as you go.
    This helps prevent you from racing ahead of learners
    and makes it easier for you to improvise in response to their questions.

*   Point learners at the online versions of the long-form notes
    (either on your bootcamp's home page
    or at [http://software-carpentry.org/v5/](http://software-carpentry.org/v5/)
    *after* the lesson is done:
    if you do it before the lesson,
    they'll try to read the notes while you're trying to talk.

*   Keep the SVG's of the diagrams handy in the directory where you're doing your teaching
    so that you can include them in your notebooks by adding an `<img src="...">` element to a Markdown cell
    (or just display them in your browser).
    Most people don't ever actually do this though,
    either because they forget to
    or because they have a whiteboard or flipchart handy.

*   There are (at least) three ways to get data files to learners at the start of a lesson:
    1.  Create a zip file, add it to your bootcamp's repository, and put a link to it in your bootcamp's `index.html` page
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
