---
layout: lesson
root: ../..
title: Introducing Version Control
level: intermediate
---
Barbara Biologist and Candace Cosmologist have been hired by Universal Missions
(a space services spinoff from Euphoric State University)
to figure out where the company should send its next planetary lander.
Since their couple biocosmology simulation software is still being developed,
they need to coordinate the changes between
Barbara's biology code and Candace's new cosmology simulations.
Barbara and Candace are developing more than just a report;
their reproducible analysis will include data, code, figures, and text.
New online writing collaboration tools
such as writeLaTeX , Google Docs, and Microsoft SkyDrive
can support their distributed authoring process,
but so far they've relied on email and Dropbox
to coordinate the rest of their development.

Although they want to be able to work on the report at the same time,
they have run into problems doing this in the past.
If they take turns,
each one will spend a lot of time waiting for the other to finish,
but if they work on their own copies and email changes back and forth,
things will be lost, overwritten, or duplicated.
Barbara and Candace also fear
this process of communicating changes
to code, data, figures, and text
won't scale as they add more scientists to their team.
A final concern,
is that the text and graphics in their report,
might fall out of sync with the code that generated them,
and that they might accidentally publish a chart,
their code can't reproduce
and that they can no longer explain.

Barbara and Candace decide to use [ADDTHISDESCRIPTION][distributed version control](../gloss.html#version-control)
to manage their workflow.
Distributed version control is better than mailing files back and forth or cloud-based syncing services because:

*   Nothing that is committed to version control is ever lost.
    This means it can be used like the "undo" feature in an editor,
    and since all old versions of files are saved
    it's always possible to go back in time to see exactly who wrote what on a particular day,
    or what version of a program was used to generate a particular set of results.
*   It keeps a record of who made what changes when,
    so that if people have questions later on,
    they know who to ask.
*   It's hard (but not impossible) to accidentally overlook or overwrite someone's changes,
    because the version control system highlights them automatically.
*   Authors can checkpoint their work, even when they're not connected to the Internet.
    This means that Barbara can experimentally try out new predator-prey models
    without worrying about changing the main simulation,
    and that Candace can checkpoint her work on the report at any time,
    even when her cable provider's Internet connection,
    mysteriously fails in the middle of the day.

This lesson shows how to use
a popular open source distributed version control system called Git.
It is more complex than some alternatives,
but it is widely used,
and is further augmented by a free hosting site called [GitHub](http://github.com).
No matter which distributed version control system you use,
the most important thing to learn is not the details of their more obscure commands,
but the fundamentals of how they work,
and the workflow that they encourage.
