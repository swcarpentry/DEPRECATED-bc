---
layout: lesson
root: ../..
title: Automating an analysis pipeline using doit
level: intermediate
---

This lesson series covers the automation of analysis pipeline's using
python's doit library. You should be familiar with functions, libraries
and dictionaries in python. Additionally, knowledge of python generators
would be very helpful.

Complex data analysis often involves a series of steps that all have 
to be carried out in a specific order, and may also require the
creation of a number of intermediate files. These types of pipelines
are well suited to automation using a class of software tools called
"build tools". Although this lesson covers one specific tool, many
of the key concepts should apply to other build tools. If you aren't 
a big python fan, there will almost certainly be tools available in
the language of your choice.

The most popular build tool is Make. For a specific comparison between
Make and doit, see [Make vs. doit](make-vs-doit.html).

<div class="toc" markdown="1">

1.  [Getting started with doit](01-doit_basics.html)
2.  [Sub-tasks](02-sub_tasks.html)
3.  [Checking whether analyses are up to date](03-uptodate.html)

</div>
