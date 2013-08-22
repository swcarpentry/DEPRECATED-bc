---
layout: lesson
root: ../..
title: Instructor's Guide for Version Control with Git
---

## Opening

Wolfman and Dracula have been hired by Universal Missions
(a space services spinoff from Euphoric State University)
to figure out where the company should send its next planetary lander.
They want to be able to work on the plans at the same time,
but they have run into problems doing this in the past.
If they take turns,
each one will spend a lot of time waiting for the other to finish.
On the other hand,
if they work on their own copies and email changes back and forth
they know that things will be lost, overwritten, or duplicated.

The right solution is to use [version control](glossary.html#version_control)
to manage their work.
Version control is better than mailing files back and forth because:

- It's hard (but not impossible) to accidentally overlook or overwrite someone's changes,
  because the version control system highlights them automatically.
- It keeps a record of who made what changes when,
  so that if people have questions later on,
  they know who to ask
  (or blame).
- Nothing that is committed to version control is ever lost.
  This means it can be used like the "undo" feature in an editor,
  and since all old versions of files are saved
  it's always possible to go back in time to see exactly who wrote what on a particular day,
  or what version of a program was used to generate a particular set of results.

The rest of this chapter will explore how to use
a popular open source version control system called Git.
It is more complex than some older systems like Subversion,
but it is widely used,
primarily because of a hosting site called [GitHub](http://github.com).
No matter which system you use,
the most important thing to learn is not the details of their more obscure commands,
but the workflow that they encourage.

## For Instructors

Version control is the most important practical skill we introduce.
Provided there aren't network problems,
each half of the lesson can be covered in about an hour,
but don't rush:
as the last paragraph of the introduction above says,
the workflow matters more than the ins and outs of any particular tool.
By the end,
the instructor should be able to get learners to chant,
"Branch, edit, commit, merge" in unison,
and have them understand what those terms mean
and why that's a good way to structure their working day.

Here's some guidance for teaching this material:

- Make sure the network is working <em>before</em> starting this lesson.
- Give learners a ten-minute overview of what version control does for them
  before diving into the watch-and-do practicals.
  Most of them will have tried to co-author papers by emailing files back and forth,
  or will have biked into the office
  only to realize that the USB key with last night's work
  is still on the kitchen table.
  Instructors can also make jokes about directories with names like
  "final version",
  "final version revised",
  "final version with reviewer three's corrections",
  "really final version",
  and,
  "come on this really has to be the last version"
  to motivate version control as a better way to collaborate
  and as a better way to back work up.
- Version control is typically taught after the shell,
  so collect learners' names during that session
  and create a repository for them to share
  with their names as both their IDs and their passwords.
  The easiest way to create the repository is to use GitHub;
  if they don't already understand SSH,
  stick to using the HTTP protocol for everything.
- Be very clear what files learners are to edit
  and what user IDs they are to use
  when giving instructions.
  It is common for them to edit the instructor's biography,
  or to use the instructor's user ID and password when committing.
  Be equally clear *when* they are to edit things:
  it's also common for someone to edit the file the instructor is editing
  and commit changes while the instructor is explaining what's going on,
  so that a conflict occurs when the instructor comes to commit the file.
- Learners could do most exercises with repositories on their own machines,
  but it's hard for them to see how version control helps collaboration
  unless they're sharing a repository with other learners.

## Prerequisites

Basic shell concepts and skills (`ls`, `cd`, `mkdir`, editing files).

## A Better Backup

Git was built by programmers for programmers,
which means that like many such tools
it has a bewildering variety of knobs and dials.
To get started with it,
let's open a shell and configure a few things:

```
$ git config --global user.name "Vlad Dracula"
$ git config --global user.email "vlad@tran.sylvan.ia"
$ git config --global color.ui "auto"
```

Git commands are written `git verb`,
where `verb` is what we actually want it to do.
In this case,
we're setting three global configuration variables to tell it
our name,
our email address,
and that we want it to colorize output.

We can now start actually using Git.
Let's create a directory for our work:

```
$ mkdir planets
$ cd planets
```

and tell Git to initialize it:

```
$ git init .
```

If we use `ls` to show the directory's contents,
it appears that nothing has changed:

```
$ ls
```

But if we add the `-a` flag to show everything,
we can see that Git has created a hidden directory called `.git`:

```
$ ls -a
.	..	.git
```

Git will store information about our project in this directory.
If you ever delete it,
you will lose the history of your project,
so please don't.

We can ask Git for the status of our project at any time like this:

```
$ git status
# On branch master
#
# Initial commit
#
nothing to commit (create/copy files and use "git add" to track)
```

Let's add some notes about Mars's suitability as a base.
(We'll echo the text to the file so that you can see what we're doing,
but in real life you would use a text editor.)

```
$ echo "Cold and dry, but everything is my favorite color" > mars.txt
$ ls
mars.txt
$ git status
# On branch master
#
# Initial commit
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	mars.txt
nothing added to commit but untracked files present (use "git add" to track)
```

The message "untracked files" means that there's a file in the directory
that Git doesn't think it's repsonsible for managing.
We can tell it that it should start like this:

```
$ git add mars.txt
```

and check that the right thing happened like this:

```
$ git status
# On branch master
#
# Initial commit
#
# Changes to be committed:
#   (use "git rm --cached <file>..." to unstage)
#
#	new file:   mars.txt
#
```

Git now knows that it's supposed to keep tack of this file,
but it *hasn't* recorded our changes for posterity---not yet.
To do that,
we need to run one more command:

```
$ git commit -m "Starting to think about Mars"
[master (root-commit) f22b25e] Starting to think about Mars
 1 file changed, 1 insertion(+)
 create mode 100644 mars.txt
```

When we run `git commit`,
Git takes everything we have told it to save
and stores a copy permanently inside its special `.git` directory
so that we can recover it later if we want to.
We use the `-m` flag to specify a comment that we want saved as well
to help us remember later on what we did and why.
We can use `git status` to check that everything has been saved:

```
$ git status
# On branch master
nothing to commit, working directory clean
```

We'll come back and explain what `branch master` means soon;
for the moment,
all we need to know is that once Git has saved things,
we can ask it about their history:

```
$ git log
commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 09:51:46 2013 -0400

    Starting to think about Mars
```

Now suppose Dracula adds more information to the file
(remember, `>>` appends rather than overwriting):

```
$ echo "The two moons may be a problem for Wolfman" >> mars.txt
```

This time, `git status` tells us that the file has been modified,
because Git already knows it's supposed to keep track of it:

```
$ git status
# On branch master
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
no changes added to commit (use "git add" and/or "git commit -a")
```

The key phrase is in the last line:
"no changes added to commit".
We have changed this file,
but we haven't committed to making those changes yet.
Let's double-check our work using `git diff`,
which shows us the differences between
the current state of the file
and the most recently saved version:

```
$ git diff
diff --git a/mars.txt b/mars.txt
index df0654a..315bf3a 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,2 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
```

The output is rather cryptic,
but we can break it down into pieces:

1. The first line tells us that Git is using the Unix `diff` command
   to compare the old and new versions of the file.
2. The second line tells exactly which versions of the file it is comparing;
   we'll look in a moment at what `df0654a` and `315bf3a` mean.
3. The remaining lines show us the actual differences
   and the lines on which they occur.
   The numbers between the `@@` markers tell editors which lines we're changing,
   and if you look in the left margin below them,
   you'll see the line we are adding marked with a '+'.

Let's commit our change:

```
$ git commit -m "Concerns about Mars's moons on my furry friend"
# On branch master
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
no changes added to commit (use "git add" and/or "git commit -a")
```

Whoops:
Git refuses to commit the changes because we didn't use `git add` first.
Let's do that:

```
$ git add mars.txt
$ git commit -m "Concerns about Mars's moons on my furry friend"
[master 34961b1] Concerns about Mars's moons on my furry friend
 1 file changed, 1 insertion(+)
```

Git insists that we add files to the set we want to commit
before actually committing anything
because we frequently won't want to commit everything at once.
For example,
suppose we're adding a few more citations to our supervisor's work
to our thesis.
We might want to commit those additions,
and the corresponding addition to the bibliography,
but *not* commit the work we've been doing on the conclusion.
To allow for this,
Git has a special staging area
where it keeps track of things that have been added to
the current [change set](glossary.html#change_set)
but not yet committed.
`git add` puts things in this area,
and `git commit` then copies them to long-term storage:

FIXME: diagram

The following commands show this in action:

```
$ echo "But the Mummy will appreciate the lack of humidity" >> mars.txt
$ git diff
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
```

So far, so good:
we've made a change,
and `git diff` tells us what it is.
Now let's put that change in the staging area
and see what `git diff` reports:

```
$ git add mars.txt
$ git diff
```

There is no output:
as far as Git can tell,
there's no difference between what it's been asked to record
and what's currently in the directory.
However,
if we do this:

```
$ git diff --staged
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
```

it shows us the difference between
the last committed change
and what's in the staging area.
Let's save our changes:

```
$ git commit -m "Thoughts about the climate"
[master 005937f] Thoughts about the climate
 1 file changed, 1 insertion(+)
```

check our status:

```
$ git status
# On branch master
nothing to commit, working directory clean
```

and look at the history of what we've done so far:

```
$ git log
git log
commit 005937fbe2a98fb83f0ade869025dc2636b4dad5
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 10:14:07 2013 -0400

    Thoughts about the climate

commit 34961b159c27df3b475cfe4415d94a6d1fcd064d
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 10:07:21 2013 -0400

    Concerns about Mars's moons on my furry friend

commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 09:51:46 2013 -0400

    Starting to think about Mars
```

If we want to see what we changed when,
we can use `git diff` yet again.
We can refer to old versions
using the notation `HEAD~1`, `HEAD~2`, and so on:

```
$ git diff HEAD~1 mars.txt
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity

$ git diff HEAD~2 mars.txt
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
```

`HEAD` is a special term meaning "the most recently saved version".
`HEAD~1` is pronounced "head minus one",
and means "the previous revision".
We can also refer to revisions using
those long strings of digits and letters
that `git log` displays.
These are unique IDs for the changes,
and "unique" really does mean unique:
every change to any set of files on any machine
has a unique 40-character identifier.
Our first commit was given the ID
f22b25e3233b4645dabd0d81e651fe074bd8e73b,
so let's try this:

```
$ git diff f22b25e3233b4645dabd0d81e651fe074bd8e73b mars.txt
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
```

That's the right answer,
but typing in 40-character strings is annoying,
so Git lets us use just the first few:

```
$ git diff f22b25e mars.txt
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
```

```
echo "### Show how to undo recent changes."
echo "Fourth line of first file" >> first.txt
cat first.txt
git reset --hard HEAD
cat first.txt
git log
```

```
echo "### Show how to undo deeper changes."
cat first.txt
git reset --hard HEAD~1
cat first.txt
git log
```

## Branching

```
echo "### Set up repository."
mkdir alpha
cd alpha
git init
echo "First line" > first.txt
git add first.txt
git commit -a -m "Setting up."

echo "### Creating a branch"
git branch analysis
git branch

echo "### Making changes on main branch."
echo "Second line added in master" >> first.txt
git add first.txt
git commit -m "Making a change in master"
git log
ls

echo "### Making changes on analysis branch."
git checkout analysis
ls
git branch
echo "Second file in analysis" > second.txt
git add second.txt
git commit -m "Adding a file in analysis"
git log
ls

echo "### Bring that change into the master."
git checkout master
ls
git merge -m "Bringing in stuff from 'analysis' branch" analysis
ls
git status
```

```
git log --oneline --topo-order --graph
```

## Cloning

```
echo "### Set up master copy."
mkdir alpha
pushd alpha
git init --bare
popd
```

```
echo "### Create clones."
git clone alpha beta
git clone alpha gamma
```

```
echo "### Push some content."
pushd beta
git config push.default simple
git status
git remote -v
echo "First line" > first.txt
git add first.txt
git commit -m "Adding first file"
git status
git push origin
popd
```

```
echo "### Pull content."
pushd gamma
git config push.default simple
git remote -v
git pull origin
ls
cat first.txt
popd
```

```
echo "### Show conflicts."
pushd beta
echo "Second line from beta" >> first.txt
git add -A
git commit -m "Beta's copy of first"
git push origin
popd
pushd gamma
echo "Second line from gamma" >> first.txt
git add -A
git commit -m "Gamma's copy of first"
git pull origin
echo "First line" > first.txt
echo "Second line from beta" >> first.txt
echo "Second line from gamma" >> first.txt
cat first.txt
git commit -a -m "Resolving merge conflict"
git log
git push origin
popd
pushd beta
git pull origin
cat first.txt
popd
```

## Collaborating

FIXME
