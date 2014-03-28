---
layout: lesson
root: ../..
title: Version Control with Git
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
  Most of them will have tried to co-author papers by emailing files back and forth
  (or in more recent years using cloud storage with file synchronization,
  e.g. Dropbox, Google Drive, ...),
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

    $ git config --global user.name "Vlad Dracula"
    $ git config --global user.email "vlad@tran.sylvan.ia"
    $ git config --global color.ui "auto"

Git commands are written `git verb`,
where `verb` is what we actually want it to do.
In this case,
we're setting three global configuration variables to tell it
our name,
our email address,
and that we want it to colorize output.

There are a few more configuration variables to set depending on your OS. First, 
choose a text editor. We recommend that novices use [GNU
nano](http://en.wikipedia.org/wiki/GNU_nano) because it's easy to use and
works in most operating systems. Some other options
might be TextEdit on the Mac, gedit on GNU/Linux or Notepad on Windows. The
default on many systems is vi, which is not a friendly text editor for beginners.
If they've installed a better editor for the workshop, use that instead.
Make sure the editor runs from the command line as configured; use a full path if necessary. 

    $ git config --global core.editor "nano"

The second OS-specific option deals with the different handling of line endings. If they ever collaborate 
with a computer on another OS, this configuration will prevent headaches.

On Mac:

    $ git config --global core.autocrlf "input"


On GNU/Linux:

    $ git config --global core.autocrlf "input"


On Windows:

    $ git config --global core.autocrlf "true"


We can now start actually using Git.
Let's create a directory for our work:

    $ mkdir planets
    $ cd planets

and tell Git to initialize it:

    $ git init

If we use `ls` to show the directory's contents,
it appears that nothing has changed:

    $ ls

But if we add the `-a` flag to show everything,
we can see that Git has created a hidden directory called `.git`:

    $ ls -a
    .	..	.git

Git will store information about our project in this directory.
If you ever delete it,
you will lose the history of your project,
so please don't.

We can ask Git for the status of our project at any time like this:

    $ git status
    # On branch master
    #
    # Initial commit
    #
    nothing to commit (create/copy files and use "git add" to track)

Let's add some notes about Mars's suitability as a base.
(We'll echo the text to the file so that you can see what we're doing,
but in real life you would use a text editor.)

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

The message "untracked files" means that there's a file in the directory
that Git doesn't think it's repsonsible for managing.
We can tell it that it should start like this:

    $ git add mars.txt

and check that the right thing happened like this:

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

Git now knows that it's supposed to keep track of this file,
but it *hasn't* recorded our changes for posterity---not yet.
To do that,
we need to run one more command:

    $ git commit -m "Starting to think about Mars"
    [master (root-commit) f22b25e] Starting to think about Mars
     1 file changed, 1 insertion(+)
     create mode 100644 mars.txt

When we run `git commit`,
Git takes everything we have told it to save
and stores a copy permanently inside its special `.git` directory
so that we can recover it later if we want to.
We use the `-m` flag to specify a comment that we want saved as well
to help us remember later on what we did and why.
We can use `git status` to check that everything has been saved:

    $ git status
    # On branch master
    nothing to commit, working directory clean

We'll come back and explain what `branch master` means soon;
for the moment,
all we need to know is that once Git has saved things,
we can ask it about their history:

    $ git log
    commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
    Author: Vlad Dracula <vlad@tran.sylvan.ia>
    Date:   Thu Aug 22 09:51:46 2013 -0400
    
        Starting to think about Mars

Now suppose Dracula adds more information to the file
(remember, `>>` appends rather than overwriting):

    $ echo "The two moons may be a problem for Wolfman" >> mars.txt

This time, `git status` tells us that the file has been modified,
because Git already knows it's supposed to keep track of it:

    $ git status
    # On branch master
    # Changes not staged for commit:
    #   (use "git add <file>..." to update what will be committed)
    #   (use "git checkout -- <file>..." to discard changes in working directory)
    #
    #	modified:   mars.txt
    #
    no changes added to commit (use "git add" and/or "git commit -a")

The key phrase is in the last line:
"no changes added to commit".
We have changed this file,
but we haven't committed to making those changes yet.
Let's double-check our work using `git diff`,
which shows us the differences between
the current state of the file
and the most recently saved version:

    $ git diff
    diff --git a/mars.txt b/mars.txt
    index df0654a..315bf3a 100644
    --- a/mars.txt
    +++ b/mars.txt
    @@ -1 +1,2 @@
     Cold and dry, but everything is my favorite color
    +The two moons may be a problem for Wolfman

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

    $ git commit -m "Concerns about Mars's moons on my furry friend"
    # On branch master
    # Changes not staged for commit:
    #   (use "git add <file>..." to update what will be committed)
    #   (use "git checkout -- <file>..." to discard changes in working directory)
    #
    #	modified:   mars.txt
    #
    no changes added to commit (use "git add" and/or "git commit -a")

Whoops:
Git refuses to commit the changes because we didn't use `git add` first.
Let's do that:

    $ git add mars.txt
    $ git commit -m "Concerns about Mars's moons on my furry friend"
    [master 34961b1] Concerns about Mars's moons on my furry friend
     1 file changed, 1 insertion(+)

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

So far, so good:
we've made a change,
and `git diff` tells us what it is.
Now let's put that change in the staging area
and see what `git diff` reports:

    $ git add mars.txt
    $ git diff

There is no output:
as far as Git can tell,
there's no difference between what it's been asked to record
and what's currently in the directory.
However,
if we do this:

    $ git diff --staged
    diff --git a/mars.txt b/mars.txt
    index 315bf3a..b36abfd 100644
    --- a/mars.txt
    +++ b/mars.txt
    @@ -1,2 +1,3 @@
     Cold and dry, but everything is my favorite color
     The two moons may be a problem for Wolfman
    +But the Mummy will appreciate the lack of humidity

it shows us the difference between
the last committed change
and what's in the staging area.
Let's save our changes:

    $ git commit -m "Thoughts about the climate"
    [master 005937f] Thoughts about the climate
     1 file changed, 1 insertion(+)

check our status:

    $ git status
    # On branch master
    nothing to commit, working directory clean

and look at the history of what we've done so far:

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

If we want to see what we changed when,
we can use `git diff` yet again.
We can refer to old versions
using the notation `HEAD~1`, `HEAD~2`, and so on:

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

    $ git diff f22b25e3233b4645dabd0d81e651fe074bd8e73b mars.txt
    diff --git a/mars.txt b/mars.txt
    index df0654a..b36abfd 100644
    --- a/mars.txt
    +++ b/mars.txt
    @@ -1 +1,3 @@
     Cold and dry, but everything is my favorite color
    +The two moons may be a problem for Wolfman
    +But the Mummy will appreciate the lack of humidity

That's the right answer,
but typing in 40-character strings is annoying,
so Git lets us use just the first few:

    $ git diff f22b25e mars.txt
    diff --git a/mars.txt b/mars.txt
    index df0654a..b36abfd 100644
    --- a/mars.txt
    +++ b/mars.txt
    @@ -1 +1,3 @@
     Cold and dry, but everything is my favorite color
    +The two moons may be a problem for Wolfman
    +But the Mummy will appreciate the lack of humidity

All right:
we can save changes to files and see what we've changed---how
can we restore older versions of things?
Let's suppose we accidentally overwrite our file
by using `>` instead of `>>`:

    $ echo "We will need to manufacture our own oxygen" > mars.txt
    $ cat mars.txt
    We will need to manufacture our own oxygen

`git status` now tells us that the file has been changed,
but those changes haven't been staged:

    $ git status
    # On branch master
    # Changes not staged for commit:
    #   (use "git add <file>..." to update what will be committed)
    #   (use "git checkout -- <file>..." to discard changes in working directory)
    #
    #	modified:   mars.txt
    #
    no changes added to commit (use "git add" and/or "git commit -a")

We can put things back the way they were like this:

    $ git reset --hard HEAD
    HEAD is now at 005937f Thoughts about the climate
    $ cat mars.txt
    Cold and dry, but everything is my favorite color
    The two moons may be a problem for Wolfman
    But the Mummy will appreciate the lack of humidity

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

FIXME: diagram

Let's run this command:

    $ git branch moons

It appears to do nothing,
but behind the scenes,
it has created a new [branch](glossary.html#branch) called `moons`:

    $ git branch
    * master
      moons

FIXME: diagram

Git is now maintaining two named bookmarks in our history:
`master`,
which was created automatically when we set up the repository,
and `moons`,
which we just made.
They both point to the same revision right now,
but we can change that.
Let's make `moons` the active branch:

    $ git checkout moons
    Switched to branch 'moons'
    $ git branch
      master
    * moons

Our file looks the same:

    $ cat mars.txt
    Cold and dry, but everything is my favorite color
    The two moons may be a problem for Wolfman
    But the Mummy will appreciate the lack of humidity

because it *is* the same:
Git hasn't made a copy of it yet
because it hasn't needed to.
Let's add another line to it:

    $ echo "Maybe we should put the base on one of the moons instead?" >> mars.txt

and add an entirely new file:

    $ echo "Phobos is larger than Deimos" > moons.txt
    $ ls
    mars.txt    moons.txt

Git now tells us that we have one changed file and one new file:

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

Let's add and commit those changes
(the `-A` flag to `git add` means "add everything"):

    $ git add -A .
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

Our repository is now in the state shown below:

FIXME: diagram

The `moons` branch has advanced to record the changes we just made,
but `master` is still where it was.
If we switch back to `master`:

    $ git checkout master

our changes seem to disappear:

    $ ls
    mars.txt
    $ cat mars.txt
    Cold and dry, but everything is my favorite color
    The two moons may be a problem for Wolfman
    But the Mummy will appreciate the lack of humidity

They're still in the repository---they're just not in
the revision that `master` is currently pointing to.
In essence,
we've created a parallel timeline
that shares some history with the original one
before diverging.

Let's make some changes in the `master` branch
to further illustrate this point:

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

Our repository is now in the state shown below:

FIXME: diagram

Both `master` and `moons` have moved on from their original common state.
They could continue independent existence indefinitely,
but at some point we'll probably want to [merge](glossary.html#merge) our changes.
Let's do that now:

    $ git branch
    * master
      moons
    $ git merge moons

When we run the `git merge` command,
Git opens an editor to let us write a log entry about what we're doing.
The editor session initially contains this:

    Merge branch 'moons'
    
    # Please enter a commit message to explain why this merge is necessary,
    # especially if it merges an updated upstream into a topic branch.
    #
    # Lines starting with '#' will be ignored, and an empty message aborts
    # the commit.

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

    Merge made by the 'recursive' strategy.
     mars.txt  | 1 +
     moons.txt | 1 +
     2 files changed, 2 insertions(+)
     create mode 100644 moons.txt

We now have all of our changes in one place:

    $ ls
    mars.txt    moons.txt    names.txt

and our repository looks like this:

FIXME: diagram

We can visualize the history with this set of arguments to `git log`:

    $ git log --oneline --topo-order --graph
    *   e0cf8ab Merge branch 'moons'
    |\  
    | * 62e7791 Thinking about the moons
    * | dfcf908 We will need a cool name for our secret base
    |/  
    * 005937f Thoughts about the climate
    * 34961b1 Concerns about Mars's moons on my furry friend
    * f22b25e Starting to think about Mars

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

FIXME: diagram

Now let's switch to the `moons` branch and make a different change there:

    $ git checkout moons
    $ echo "This line added in the moons branch" >> moons.txt
    $ cat moons.txt
    Phobos is larger than Deimos
    This line added in the moons branch
    $ git add moons.txt
    $ git commit -m "Adding a line in the moons branch"
    [moons 07ebc69] Adding a line in the moons branch
     1 file changed, 1 insertion(+)

Our repository now looks like this:

FIXME: diagram

Let's pull all the changes made in `master` into the `moons` branch:

    $ git merge master
    Auto-merging moons.txt
    CONFLICT (content): Merge conflict in moons.txt
    Automatic merge failed; fix conflicts and then commit the result.

Git has detected that the changes made in the `master` branch
overlap with those made in the `moons` branch.
If we ask it for our status,
we get this:

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

which tells us that it brought over the file `names.txt` successfully
(which was added in `master`, but didn't yet exist in `moons`),
but was unable to handle the conflict in `moons.txt`.
What it *has* done is mark the conflict in that file:

    $ cat moons.txt
    Phobos is larger than Deimos
    <<<<<<< HEAD
    This line added in the moons branch
    =======
    This line added in master
    >>>>>>> master

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

    $ cat moons.txt
    Phobos is larger than Deimos
    Lines added in the master and moons branches

To finish merging,
we need to add `moons.txt` to the changes being made by the merge
and then commit:

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

Our repository now looks like this:

FIXME: diagram

Git tries hard to keep track of what we've merged with what,
so if we switch back to `master` and merge the changes in `moons`,
we don't have to fix things by hand again:

    $ git checkout master
    $ git merge moons
    Updating 5ae9631..2f20801
    Fast-forward
     moons.txt | 2 +-
     1 file changed, 1 insertion(+), 1 deletion(-)
    $ cat moons.txt 
    Phobos is larger than Deimos
    Lines added in the master and moons branches

The key phrase here is "fast-forward"
(which appears in the output of the `git merge` command).
Since we had already resolved the conflicts between
the copies of `moons.txt` in the `master` and `moons` branches,
Git brings the result over on its own.

## Collaborating

Version control really comes into its own
when we begin to collaborate with other people.
We already have most of the machinery we need to do this:
repositories,
branches,
and the `commit` and `merge` commands.
The last trick is to merge from branches that are in other repositories,
not our own.

Systems like Git and Mercurial allow us to merge changes
between any two repositories.
In practice,
though,
it's easiest to use a definitive master copy as a central hub,
and for that master copy to be on the web rather than on someone's laptop
(so that it's accessible even when that "someone" is off the network).
Most programmers use hosting services like [GitHub](http://github.com) or [BitBucket](http://bitbucket.org)
to hold those master copies;
we'll explore the pros and cons of this in the final section of this lesson,
but will use GitHub until then.

Let's start by sharing the changes we've made to our current project with the world.
Log in to GitHub,
then create a new repository called `planets`
using their GUI:

FIXME: screenshot

This effectively does the following on GitHub's servers:

    $ mkdir planets
    $ cd planets
    $ git init

We're now in the situation shown in the figure below:

FIXME: diagram

Our local repository still has two branches called `master` and `moons`,
with the same contents as before.
The remote repository on GitHub only has a single branch,
`master`,
and doesn't contain any files yet.

The next step---the crucial one---is to connect the two repositories.
We do this by making the GitHub repository a [remote](glossary.html#remote_repository)
for the local repository.
The home page of the repository on GitHub includes
the string we need to identify it:

FIXME: screenshot

For now,
we'll use the 'http' identifier,
since it requires the least setup.
Copy that string from the browser,
go into the local `planets` repository,
and run this command:

    $ git remote add origin https://github.com/yourname/planets

(using your GitHub ID instead of `yourname`).
We can check that the command has worked by running `git remote -v`:

    $ git remote -v
    origin   https://github.com/yourname/planets.git (push)
    origin   https://github.com/yourname/planets.git (fetch)

There's nothing special about the name `origin`:
we can use almost anything,
but we'll see in a moment why `origin` is a sensible choice.
Once this is set up,
the following command will push the changes from our local repository's `master` branch
to the corresponding branch in the repository on GitHub:

    $ git push origin master
    Counting objects: 27, done.
    Delta compression using up to 4 threads.
    Compressing objects: 100% (23/23), done.
    Writing objects: 100% (27/27), 2.62 KiB, done.
    Total 27 (delta 5), reused 0 (delta 0)
    To https://github.com/gvwilson/planets.git
     * [new branch]      master -> master

This command just did what `git merge` does,
except it moved changes between repositories
rather than just between branches.
Our local and remote repositories are now in this state:

FIXME: diagram

We can pull changes from the remote repository to the local one as well:

    $ git pull origin master
    From https://github.com/gvwilson/planets
     * branch            master     -> FETCH_HEAD
    Already up-to-date.

Pulling has no effect in this case
because the two repositories are already synchronized.
If someone else had pushed some changes,
though,
this command would download them to our local repository:

FIXME: diagram

The model shown above,
in which everyone pushes and pulls from a single repository,
is perfectly usable,
but there's one thing it *doesn't* let us do,
and that's [code review](glossary.html#code_review).
Suppose Dracula wants to be able to look at Wolfman's changes
before merging them into the master copy on GitHub,
just as he would review Wolfman's paper before publishing it
(or perhaps even before submitting it for publication).
We need to arrange things so that Wolfman can commit his changes
and Dracula can compare them with the master copy;
in fact,
we want Wolfman to be able to commit many times,
so that he can incorporate Dracula's feedback
and get further review
as often as necessary.

Rather than the model shown above,
most programmers therefore use a slightly more complex model.
When the project starts,
Dracula creates a repository on GitHub
in exactly the same way as we created the `planets` repository a few moments ago,
and then [clones](glossary.html#repository_clone) it to his desktop:

    $ git clone https://github.com/vlad/undersea.git
    Cloning into 'undersea'...
    warning: You appear to have cloned an empty repository.

`git clone` automatically adds the original repository on GitHub
as a remote of the local repository called `origin`---this
is why we chose `origin` as a remote name in our previous example:

    $ cd undersea
    $ git remote -v
    origin	    https://github.com/vlad/undersea.git (fetch)
    origin	    https://github.com/vlad/undersea.git (push)

Dracula can now push and pull changes just as before.

Wolfman doesn't clone Dracula's GitHub repository directly.
Instead,
he [forks](glossary.html#fork_repository) it,
i.e.,
clones it on GitHub.
He does this using the GitHub web interface:

FIXME: screenshot

He then clones his GitHub repository,
not Dracula's,
to give himself a desktop copy:

FIXME: diagram

This may seem like unnecessary work,
but it allows Wolfman and Dracula to collaborate much more effectively.
Suppose Wolfman makes a change to the project.
He commits it to his local repository,
then uses `git push` to copy those changes to GitHub:

FIXME: diagram

He then creates a [pull request](glossary.html#pull_request),
which notifies Dracula that Wolfman wants to merge some changes
into Dracula's repository:

FIXME: screenshot

A pull request is a merge waiting to happen.
When Dracula views it online,
he can see and comment on the changes Wolfman wants to make:

FIXME: screenshot

Commenting is the crucial step here,
and half the reason Wolfman went to the trouble of forking the repository on GitHub.
Dracula,
or anyone else involved in the project,
can now give Wolfman feedback on what he is trying to do:
this function is too long,
that one contains a bug,
there's a special case that isn't being handled anywhere,
and so on.
Wolfman can then update his code,
commit locally,
and push those changes to GitHub
to update the pull request:

FIXME: diagram

This process is exactly like peer review of papers,
though usually much faster.
In large open source projects like Firefox,
it's very common for a pull request to be updated several times
before finally being accepted and merged.
Working this way not only helps maintain the quality of the code,
it is also a very effective way to transfer knowledge.

What happens if Wolfman wants to do more work
while he's waiting for Dracula to review his first modification?
Simple:
he creates a new branch in his local repository,
pushes it to GitHub,
and then issues a pull request from that:

FIXME: diagram

We can now see why Git, Mercurial, and other modern version control systems
use branching so much:
it helps people work concurrently but asynchronously,
i.e.,
in parallel but on their own time.
It might take Dracula several days to get around to reviewing Wolfman's changes.
Rather than being stalled until then,
Wolfman can just switch to another branch and work on something else,
then switch back when Dracula's review finally comes in.
Once the changes in a particular branch have been accepted,
Wolfman can delete it;
provided it has been merged into `master` (or some other branch),
the only thing that will be lost is the pointer with the branch name,
not the changes themselves.

We said above that code review is
half the reason every developer should have their own repository on GitHub
(or whatever service is being used).
The other reason is that working this way allows people to explore ideas
without needing permission from any central authority.
If you want to change this tutorial,
you can fork the Software Carpentry repository on GitHub
and start rewriting things.
If you think we might prefer your version to ours,
you can send us a pull request,
but you don't have to.
If other people like your version better than ours,
they can start forking your repository
and sending pull requests to you instead of to us.

If this sounds familiar,
it's because it is the way science itself works.
When someone publishes a new method or result,
other scientists can immediately start building on top of it---essentially,
they can create their own fork of the work
and start committing changes to it.
If the first scientist likes the second's work,
she can incorporate those findings into her next paper,
which is analogous to merging a pull request.
If she doesn't,
then it's up to other scientists to decide whose work to build on,
or whether to try to combine both approaches.

## The Opposite of "Open" Isn't "Closed", It's "Broken"

Free sharing of information might be the ideal in science,
but the reality is often more complicated.
Normal practice today looks something like this:

* A scientist collects some data and stores it on a machine
  that is occasionally backed up by her department.
* She then writes or modifies a few small programs
  (which also reside on her machine)
  to analyze that data.
* Once she has some results,
  she writes them up and submits her paper.
  She might include her data---a growing number of journals require this---but
  she probably doesn't include her code.
* Time passes.
* The journal sends her reviews written anonymously by a handful of other people in her field.
  She revises her paper to satisfy them,
  during which time she might also modify the scripts she wrote earlier,
  and resubmits.
* More time passes.
* The paper is eventually published.
  It might include a link to an online copy of her data,
  but the paper itself will be behind a paywall:
  only people who have personal or institutional access
  will be able to read it.

For a growing number of scientists,
though,
the process looks like this:

* The data that the scientist collects is stored in an open access repository
  like [figshare](http://figshare.com/) or [Dryad](http://datadryad.org/)
  as soon as it's collected,
  and given its own DOI.
* The scientist creates a new repository on GitHub to hold her work.
* As she does her analysis,
  she pushes changes to her scripts
  (and possibly some output files)
  to that repository.
  She also uses the repository for her paper;
  that repository is then the hub for collaboration with her colleagues.
* When she's happy with the state of her paper,
  she posts a version to [arXiv](http://arxiv.org/)
  or some other preprint server
  to invite feedback from peers.
* Based on that feedback,
  she may post several revisions
  before finally submitting her paper to a journal.
* The published paper includes links to her preprint
  and to her code and data repositories,
  which  makes it much easier for other scientists
  to use her work as starting point for their own research.

<!-- SJE: Can you give a reference to support "Studies have shown"? -->
Studies have shown that the more open model accelerates discovery,
and that more open work is more widely cited.
However,
people who want to work this way need to make some decisions
about what exactly "open" means in practice.

### Licensing

The first question is licensing.
Broadly speaking,
there are two kinds of open license for software,
and half a dozen for data and publications.
For software,
people can choose between the [GNU Public License](http://opensource.org/licenses/GPL-3.0) (GPL) on the one hand,
and licenses like the [MIT](http://opensource.org/licenses/MIT)
and [BSD](http://opensource.org/licenses/BSD-2-Clause) licenses on the other.
All of these licenses allow unrestricted sharing and modification of programs,
but the GPL is [infective](glossary.html#infective_license):
anyone who distributes a modified version of the code
(or anything that includes GPL'd code)
must make *their* code freely available as well.

Proponents of the GPL argue that this requirement is needed
to ensure that people who are benefiting from freely-available code
are also contributing back to the community.
Opponents counter that many open source projects have had long and successful lives
without this condition,
and that the GPL makes it more difficult to combine code from different sources.
At the end of the day,
what matters most is that:

1. every project have a file in its home directory
   called something like `LICENSE` or `LICENSE.txt`
   that clearly states what the license is, and
2. people use existing licenses rather than writing new ones.

The second of these points is as important as the first:
most scientists are not lawyers,
so wording that may seem sensible to a layperson
may have unintended gaps or consequences.
The [Open Source Initiative](http://opensource.org/)
maintains a list of open source licenses,
and [tl;drLegal](http://www.tldrlegal.com/) explains most of them in plain English.

When it comes to data, publications, and other "static" material,
scientists have many more options to choose from.
The good news is that an organization called [Creative Commons(http://creativecommons.org/)
has prepared a set of licenses using combinations of four basic restrictions:

* Attribution: derived works must give the original author credit for their work.
* No Derivatives: people may copy the work, but must pass it along unchanged.
* Share Alike: derivative works must license their work under the same terms as the original.
* Noncommercial: free use is allowed, but commercial use is not.

These four restrictions are abbreviated "BY", "ND", "SA", and "NC" respectively,
so "CC-BY-ND" means,
"People can re-use the work both for free and commercially,
but cannot make changes and must cite the original."
These [short descriptions](http://creativecommons.org/licenses/)
summarize the six CC licenses in plain language,
and include links to their full legal formulations.

There is one other important license that doesn't fit into this categorization.
Scientists (and other people) can choose to put material in the public domain,
which is often abbreviated "PD".
In this case,
anyone can do anything they want with it,
without needing to cite the original
or restrict further re-use.
The table below shows how the six Creative Commons licenses and PD relate to one another:

<table border="1">
  <tr>
    <td></td>
    <td colspan="7" align="center">Licenses that can be used for derivative work or adaptation</td>
  </tr>
  <tr>
    <td>Original work</td> <td>by</td> <td>by-nc</td> <td>by-nc-nd</td> <td>by-nc-sa</td> <td>by-nd</td> <td>by-sa</td> <td>pd</td>
  </tr>
  <tr>
    <td>by</td>       <td>X</td> <td>X</td> <td>X</td> <td>X</td> <td>X</td> <td>X</td> <td> </td>
  </tr>
  <tr>
    <td>by-nc</td>    <td> </td> <td>X</td> <td>X</td> <td>X</td> <td> </td> <td> </td> <td> </td>
  </tr>
  <tr>
    <td>by-nc-nd</td> <td> </td> <td> </td> <td> </td> <td> </td> <td> </td> <td> </td> <td> </td>
  </tr>
  <tr>
    <td>by-nc-sa</td> <td> </td> <td> </td> <td> </td> <td>X</td> <td> </td> <td> </td> <td> </td>
  </tr>
  <tr>
    <td>by-nd</td>    <td> </td> <td> </td> <td> </td> <td> </td> <td> </td> <td> </td> <td> </td>
  </tr>
  <tr>
    <td>by-sa</td>    <td> </td> <td> </td> <td> </td> <td> </td> <td> </td> <td>X</td> <td> </td>
  </tr>
  <tr>
    <td>pd</td>       <td>X</td> <td>X</td> <td>X</td> <td>X</td> <td>X</td> <td>X</td> <td>X</td>
  </tr>
</table>

[Software Carpentry](http://software-carpentry.org/license.html)
uses CC-BY for its lessons and the MIT License for its code
in order to encourage the widest possible re-use.
Again,
the most important thing is for the `LICENSE` file in the root directory of your project
to state clearly what your license is.
You may also want to include a file called `CITATION` or `CITATION.txt`
that describes how to reference your project;
the one for Software Carpentry states:

    To reference Software Carpentry in publications, please cite both of the following:

    Greg Wilson: "Software Carpentry: Getting Scientists to Write Better
    Code by Making Them More Productive".  Computing in Science &
    Engineering, Nov-Dec 2006.

    Greg Wilson: "Software Carpentry: Lessons Learned". arXiv:1307.5448,
    July 2013.

    @article{wilson-software-carpentry-2006,
        author =  {Greg Wilson},
        title =   {Software Carpentry: Getting Scientists to Write Better Code by Making Them More Productive},
        journal = {Computing in Science \& Engineering},
        month =   {November--December},
        year =    {2006},
    }

    @online{wilson-software-carpentry-2013,
      author      = {Greg Wilson},
      title       = {Software Carpentry: Lessons Learned},
      version     = {1},
      date        = {2013-07-20},
      eprinttype  = {arxiv},
      eprint      = {1307.5448}
    }


### Hosting

The second big question for groups that want to open up their work
is where to host their code and data.
One option is for the lab, the department, or the university to provide a server,
manage accounts and backups,
and so on.
The main benefit of this is that it clarifies who owns what,
which is particularly important if any of the material is sensitive
(i.e.,
relates to experiments involving human subjects
or may be used in a patent application).
The main drawbacks are the cost of providing the service and its longevity:
a scientist who has spent ten years collecting data
would like to be sure that data will still be available ten years from now,
but that's well beyond the lifespan of most of the grants that fund academic infrastructure.

The alternative is to use a public hosting service like [GitHub](http://github.com),
[BitBucket](http://bitbucket.org),
[Google Code](http://code.google.com),
or [SourceForge](http://sourceforge.net).
All of these allow people to create repositories through a web interface,
and also provide mailing lists,
ways to keep track of who's doing what,
and so on.
They all benefit from economies of scale and network effects:
it's easier to run one large service well
than to run many smaller services to the same standard,
and it's also easier for people to collaborate if they're using the same service,
not least because it gives them fewer passwords to remember.

However,
all of these services place some constraints on people's work.
In particular,
they give users a choice:
if they're willing to share their work with others,
it will be hosted for free,
but if they want privacy,
they have to pay.
Sharing might seem like the only valid choice for science,
but many institutions may not allow researchers to do this,
either because they want to protect future patent applications
or simply because what's new is often also frightening.

### Get the Facts

Most students don't actually know who owns their work.
Does it belong to the university?
To their supervisor?
Can they take their work with them if they leave the university,
or did something in the fine print of the registration forms they signed without reading
hand that to someone else?

Most faculty and administrators don't actually know either.
Instead,
they half-remember what the rules were five or ten years ago,
but haven't kept up with changes made since.
As a result,
many research projects are open sourced by fiat:
a student or faculty member declares that the work is covered by a particular license,
puts it into a public repository,
and waits for someone to object.

This approach has worked surprisingly well so far,
but we still suggest that people ask someone from their institution's intellectual property office
to give a lunchtime talk on these issues
*before* deciding what to do and how to do it:
It's clear that tomorrow's science and scientific computing will be more open than yesterday's,
and that this will benefit everyone,
but the fewer stumbles we make along the way,
the faster we'll get there.
