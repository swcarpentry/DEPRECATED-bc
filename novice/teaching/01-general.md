---
layout: lesson
root: ../..
title: General Advice
level: novice
---
This is a placeholder for general notes about teaching.

#### Teaching Notes

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
