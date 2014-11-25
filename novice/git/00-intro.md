---
layout: lesson
root: ../..
title: Introducing Version Control
---
Wolfman and Dracula have been hired by Universal Missions
(a space services spinoff from Euphoric State University)
to investigate if it is possible to send their next planetary lander to Mars.
They want to be able to work on the plans at the same time,
but they have run into problems doing this in the past.
If they take turns,
each one will spend a lot of time waiting for the other to finish,
but if they work on their own copies and email changes back and forth
things will be lost, overwritten, or duplicated.

The right solution is to use [version control](../../gloss.html#version-control)
to manage their work.
Version control is better than mailing files back and forth because:

*   Nothing that is committed to version control is ever lost.
    This means it can be used like the "undo" feature in an editor,
    and since all old versions of files are saved
    it's always possible to go back in time to see exactly who wrote what on a particular day,
    or what version of a program was used to generate a particular set of results.
*   It keeps a record of who made what changes when,
    so that if people have questions later on,
    they know who to ask.
*   It's hard (but not impossible) to accidentally overlook or overwrite someone's changes:
    the version control system automatically notifies users
    whenever there's a conflict between one person's work and another's.

<div class="challenges" markdown="1">

#### Challenges

On Wikipedia all changes and their authors are tracked. You can go
[here](https://en.wikipedia.org/w/index.php?title=Mars&action=history)
and you will find the history of all changes done to the article about the planet
Mars. Find the last edit done last month and look at the changes made by
clicking on the "prev" link on the left of the history entry.

</div>

This lesson shows how to use
a popular open source version control system called Git.
It is more complex than some alternatives,
but it is widely used,
both because it's easy to set up
and because of a hosting site called [GitHub](http://github.com).
No matter which version control system you use,
the most important thing to learn is not the details of their more obscure commands,
but the workflow that they encourage.
