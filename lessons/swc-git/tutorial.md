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

`HEAD` means "the most recently saved version".
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

All right:
we can save changes to files and see what we've changed---how
can we restore older versions of things?
Let's suppose we accidentally overwrite our file
by using `>` instead of `>>`:

```
$ echo "We will need to manufacture our own oxygen" > mars.txt
$ cat mars.txt
We will need to manufacture our own oxygen
```

`git status` now tells us that the file has been changed,
but those changes haven't been staged:

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

We can put things back the way they were like this:

```
$ git reset --hard HEAD
HEAD is now at 005937f Thoughts about the climate
$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
```

The `--hard` argument to `git reset` tells it to throw away local changes:
without that,
`git reset` won't destroy our work.
`HEAD` tells `git reset` that we want to put things back to
the way they were recorded in the `HEAD` revision.
(Remember,
we haven't done a `git commit` with these changes yet,
so `HEAD` is still where it was.)
We can use `git reset --hard HEAD~55` and so on
to back up to earlier revisions,
`git reset --hard 34961b1` to back up to a particular revision,
and so on.

But what if we want to recover files *without* losing the work we've done since?
For example,
what if we have added some material to the conclusion of our paper that we'd like to keep,
but we want to get back an earlier version of the introduction?
To accomplish that,
we'll need to explore branching.

## Branching

Here's where we are right now:

```
$ git log
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
$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
```

FIXME: diagram

Let's run this command:

```
$ git branch moons
```

It appears to do nothing,
but behind the scenes,
it has created a new [branch](glossary.html#branch) called `moons`:

```
$ git branch
* master
  moons
```

FIXME: diagram

Git is now maintaining two named bookmarks in our history:
`master`,
which was created automatically when we set up the repository,
and `moons`,
which we just made.
They both point to the same revision right now,
but we can change that.
Let's make `moons` the active branch:

```
$ git checkout moons
Switched to branch 'moons'
$ git branch
  master
* moons
```

Our file looks the same:

```
$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
```

because it *is* the same:
Git hasn't made a copy of it yet
because it hasn't needed to.
Let's add another line to it:

```
$ echo "Maybe we should put the base on one of the moons instead?" >> mars.txt
```

and add an entirely new file:

```
$ echo "Phobos is larger than Deimos" > moons.txt
$ ls
mars.txt    moons.txt
```

Git now tells us that we have one changed file and one new file:

```
$ git status
# On branch moons
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	moons.txt
no changes added to commit (use "git add" and/or "git commit -a")
```

Let's add and commit those changes
(the `-A` flag to `git commit` means "add everything"):

```
$ git add -A
$ git status
# On branch moons
# Changes to be committed:
#   (use "git reset HEAD <file>..." to unstage)
#
#	modified:   mars.txt
#	new file:   moons.txt
#
~/planets: git commit -m "Thinking about the moons"
[moons 62e7791] Thinking about the moons
 2 files changed, 2 insertions(+)
 create mode 100644 moons.txt
```

Our repository is now in the state shown below:

FIXME: diagram

The `moons` branch has advanced to record the changes we just made,
but `master` is still where it was.
If we switch back to `master`:

```
$ git checkout master
```

our changes seem to disappear:

```
$ ls
mars.txt
$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
```

They're still in the repository---they're just not in
the revision that `master` is currently pointing to.
In essence,
we've created a parallel timeline
that shares some history with the original one
before diverging.

Let's make some changes in the `master` branch
to further illustrate this point:

```
$ echo "Should we go with a classical name like Ares Base?" > names.txt
$ git status
# On branch master
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	names.txt
nothing added to commit but untracked files present (use "git add" to track)
$ git add names.txt
$ git commit -m "We will need a cool name for our secret base"
[master dfcf908] We will need a cool name for our secret base
 1 file changed, 1 insertion(+)
 create mode 100644 names.txt
```

Our repository is now in the state shown below:

FIXME: diagram

Both `master` and `moons` have moved on from their original common state.
They could continue independent existence indefinitely,
but at some point we'll probably want to [merge](glossary.html#merge) our changes.
Let's do that now:

```
$ git branch
* master
  moons
$ git merge moons
```

When we run the `git merge` command,
Git opens an editor to let us write a log entry about what we're doing.
The editor session initially contains this:

```
Merge branch 'moons'

# Please enter a commit message to explain why this merge is necessary,
# especially if it merges an updated upstream into a topic branch.
#
# Lines starting with '#' will be ignored, and an empty message aborts
# the commit.
```

If we notice that something is wrong
and decide not to complete the merge,
we must delete everything in the file---Git interprets an empty log message to mean,
"Don't proceed."
Otherwise,
everything that isn't marked as a comment with `#` will be saved to the log.
In this case,
we'll stick with the default log message.
When we save the file and exit the editor,
Git displays this:

```
Merge made by the 'recursive' strategy.
 mars.txt  | 1 +
 moons.txt | 1 +
 2 files changed, 2 insertions(+)
 create mode 100644 moons.txt
```

We now have all of our changes in one place:

```
$ ls
mars.txt    moons.txt    names.txt
```

and our repository looks like this:

FIXME: diagram

We can visualize the history with this set of arguments to `git log`:

```
$ git log --oneline --topo-order --graph
*   e0cf8ab Merge branch 'moons'
|\  
| * 62e7791 Thinking about the moons
* | dfcf908 We will need a cool name for our secret base
|/  
* 005937f Thoughts about the climate
* 34961b1 Concerns about Mars's moons on my furry friend
* f22b25e Starting to think about Mars
```

This ASCII art is fine for small sets of changes,
but for anything significant,
it's much better to use a proper GUI that can draw graphs using lines instead of characters.

Branching strikes most newcomers as unnecessary complexity,
particularly for single-author projects.
After all,
if we need to make some changes to a project,
what do we gain by creating parallel universes?
The answer is that branching makes it easy for us
to concentrate on one thing at a time.
Suppose we are part-way through rewriting a function that calculates spatial correlations
when we realize that the task would be easier
if our file I/O routines always stored things as complex numbers.
Most people would put the spatial correlation changes aside,
change the file I/O,
then (hopefully) come back to the original task.

The problem with this is that
we have to remember what we were doing,
even if we realize halfway through rewriting file I/O
that we should also rewrite our error handling.
It's quite common to wind up with half a dozen tasks stacked on top of one another,
and quite hard to them all straight.
Branching allows us to put what we're doing in a safe place,
solve the new problem,
then resume our earlier work.

Working this way is particularly beneficial when we have a good set of unit tests.
If we pause Task A halfway through to start work on Task B,
the odds are that the tests for Task A will fail,
or that the code won't run at all because it's in pieces.
Doing things in branches avoids this.

In practice,
most developers never actually make changes directly in the `master` branch.
Instead,
they create a new branch from it for every significant change they want to make,
then merge those branches back to `master` when the work is complete.
Returning to our hypothetical example,
we would

1. create a branch called something like `better-spatial-correlation`
   for those changes:

FIXME: diagram

2. go back to master and create another branch called `file-input-produces-complex-values`
   for *those* changes:

FIXME: diagram

3. merge `file-input-produces-complex-values` into `master`:

FIXME: diagram

4. merge `master` into `better-spatial-correlation`:

FIXME: diagram

5. and then finish work on the spatial correlation function
   and merge it all back into `master`:

FIXME: diagram

And if,
partway through this process,
our supervisor asked us to change the graph-plotting routines
to conform to the university's new style guide,
we would simply switch back to `master`,
create a branch for that,
make the changes,
produce the desired graphs,
and leave the changes parked in that branch
until we were ready to merge them.

## Merging

As soon as people can work in parallel,
someone's going to step on someone else's toes.
This will even happen with a single person:
while updating the conclusion to a paper in one branch,
for example,
we may decide to make changes to the introduction,
which we have already updated in another branch.
Version control helps us manage these [conflicts](glossary.html#conflict)
by giving us tools to [merge](glossary.html#merge) overlapping changes.

To see how merging works,
we must first create a conflict.
Let's add a line to the version of `moons.txt` in the `master` branch:

```
$ git branch
* master
  moons
$ echo "This line added in master" >> moons.txt
$ cat moons.txt
Phobos is larger than Deimos
This line added in master
$ git add moons.txt
$ git commit -m "Adding a line to moons.txt in the master branch"
[master 5ae9631] Adding a line in the master branch
 1 file changed, 1 insertion(+)
```

FIXME: diagram

Now let's switch to the `moons` branch and make a different change there:

```
$ echo "This line added in the moons branch" >> moons.txt
$ cat moons.txt
Phobos is larger than Deimos
This line added in the moons branch
$ git add moons.txt
$ git commit -m "Adding a line in the moons branch"
[moons 07ebc69] Adding a line in the moons branch
 1 file changed, 1 insertion(+)
```

Our repository now looks like this:

FIXME: diagram

Let's pull all the changes made in `master` into the `moons` branch:

```
$ git merge master
Auto-merging moons.txt
CONFLICT (content): Merge conflict in moons.txt
Automatic merge failed; fix conflicts and then commit the result.
```

Git has detected that the changes made in the `master` branch
overlap with those made in the `moons` branch.
If we ask it for our status,
we get this:

```
$ git status
# On branch moons
# You have unmerged paths.
#   (fix conflicts and run "git commit")
#
# Changes to be committed:
#
#	new file:   names.txt
#
# Unmerged paths:
#   (use "git add <file>..." to mark resolution)
#
#	both modified:      moons.txt
#
```

which tells us that it brought over the file `names.txt` successfully
(which was added in `master`, but didn't yet exist in `moons`),
but was unable to handle the conflict in `moons.txt`.
What it *has* done is mark the conflict in that file:

```
$ cat moons.txt
Phobos is larger than Deimos
<<<<<<< HEAD
This line added in the moons branch
=======
This line added in master
>>>>>>> master
```

Our change---the one in `HEAD`---is preceded by `<<<<<<<`.
Git has then inserted `=======` as a separator between the conflicting changes
and marked the end of the region with `>>>>>>>`.

It is now up to us to edit this file to remove these markers
and reconcile the changes.
We can do anything we want:
keep the change in this branch,
keep the change made in the other,
write something new to replace both,
or get rid of the change entirely.
Let's replace both so that the file looks like this:

```
$ cat moons.txt
Phobos is larger than Deimos
Lines added in the master and moons branches
```

To finish merging,
we need to add `moons.txt` to the changes being made by the merge
and then commit:

```
$ git add moons.txt
$ git status
# On branch moons
# All conflicts fixed but you are still merging.
#   (use "git commit" to conclude merge)
#
# Changes to be committed:
#
#	modified:   moons.txt
#	new file:   names.txt
#
$ git commit -m "Pulling in changes from master"
[moons 2f20801] Pulling in changes from master
```

Our repository now looks like this:

FIXME: diagram

Git tries hard to keep track of what we've merged with what,
so if we switch back to `master` and merge the changes in `moons`,
we don't have to fix things by hand again:

```
$ git checkout master
$ git merge moons
Updating 5ae9631..2f20801
Fast-forward
 moons.txt | 2 +-
 1 file changed, 1 insertion(+), 1 deletion(-)
$ cat moons.txt 
Phobos is larger than Deimos
Lines added in the master and moons branches
```

The key phrase here is "fast forward"
(which appears in the output of the `git merge` command).
Since we had already resolved the conflicts between
the copies of `moons.txt` in the `master` and `moons` branches,
Git brings the result over on its own.

----------------------------------------

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

----------------------------------------

## Collaborating

FIXME
