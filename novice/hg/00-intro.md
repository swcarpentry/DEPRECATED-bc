---
layout: lesson
root: ../..
title: Introducing Version Control
---
<div class="objectives" markdown="1">

#### Objectives
* Explain when and why you should use version control

</div>

Wolfman and Dracula have been hired by Universal Missions
(a space services spin-off of Euphoric State University)
to figure out where the company should send its next planetary lander.
They want to be able to work on the plans at the same time,
but they have run into problems doing this in the past.
If they take turns,
each one will spend a lot of time waiting for the other to finish,
but if they work on their own copies and email changes back and forth
things will be lost, overwritten, or duplicated.

The solution is to use [version control](../../gloss.html#version-control)
to manage their work.
Version control is better than mailing files back and forth because:

*   Nothing that is committed to version control is ever lost.
    This means it can be used like the "undo" feature in an editor,
    and since all old versions of files are saved
    it's always possible to go back in time to see 
    exactly who wrote what on a particular day,
    or what version of a program was used to generate a particular set of results.
*   It keeps a record of who made what changes when,
    so that if people have questions later on
    they know who to ask.
*   It's hard (but not impossible) to accidentally overlook or overwrite someone's changes:
    the version control system automatically notifies users
    whenever there's a conflict between one person's work and another's.

This lesson shows how to use
a popular open source version control system called Mercurial
(also known as hg).
It is widely used,
both because it's easy to set up
and because of a hosting site called [Bitbucket](http://bitbucket.org).
No matter which version control system you use,
the most important thing to learn is not the details of their more obscure commands,
but the workflow that they encourage.
