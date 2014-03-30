---
layout: lesson
root: ../..
title: Branching in Git
---
Here's where we are right now:

<div class="in" markdown="1">
~~~
$ git log
~~~
</div>
<div class="out" markdown="1">
~~~
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
~~~
</div>

We can draw the history of the repository like this
(we'll see in a second why there's a box called `master`):

<img src="img/git-branching-01.svg" alt="Initial State of Repository" />

Let's run this command:

<div class="in" markdown="1">
~~~
$ git branch moons
~~~
</div>

It appears to do nothing,
but behind the scenes it has created a new [branch](../../gloss.html#branch) called `moons`:

<div class="in" markdown="1">
~~~
$ git branch
~~~
</div>
<div class="out" markdown="1">
~~~
* master
  moons
~~~
</div>

<img src="img/git-branching-02.svg" alt="Immediately After Creating Branch" />

Git is now maintaining two named bookmarks in our history:
`master`,
which was created automatically when we set up the repository,
and `moons`,
which we just made.
They both point to the same revision right now,
but we can change that.
Let's make `moons` the active branch:

<div class="in" markdown="1">
~~~
$ git checkout moons
~~~
</div>
<div class="out" markdown="1">
~~~
Switched to branch 'moons'
~~~
</div>
<div class="in" markdown="1">
~~~
$ git branch
~~~
</div>
<div class="out" markdown="1">
~~~
  master
* moons
~~~
</div>

<img src="img/git-branching-03.svg" alt="After Switching to Branch" />

Our file looks the same:

<div class="in" markdown="1">
~~~
$ cat mars.txt
~~~
</div>
<div class="out" markdown="1">
~~~
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~
</div>

because it *is* the same:
Let's add another line to it:

<div class="in" markdown="1">
~~~
$ echo "Maybe we should put the base on one of the moons instead?" >> mars.txt
~~~
</div>

and add an entirely new file:

<div class="in" markdown="1">
~~~
$ echo "Phobos is larger than Deimos" > moons.txt
$ ls
~~~
</div>
<div class="out" markdown="1">
~~~
mars.txt    moons.txt
~~~
</div>

Git now tells us that we have one changed file and one new file:

<div class="in" markdown="1">
~~~
$ git status
~~~
</div>
<div class="out" markdown="1">
~~~
# On branch moons
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#    modified:   mars.txt
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#    moons.txt
no changes added to commit (use "git add" and/or "git commit -a")
~~~
</div>

Let's add and commit those changes
(the `-A` flag to `git commit` means "add everything"):

<div class="in" markdown="1">
~~~
$ git add -A
$ git status
~~~
</div>
<div class="out" markdown="1">
~~~
# On branch moons
# Changes to be committed:
#   (use "git reset HEAD <file>..." to unstage)
#
#    modified:   mars.txt
#    new file:   moons.txt
#
~~~
</div>
<div class="in" markdown="1">
~~~
$ git commit -m "Thinking about the moons"
~~~
</div>
<div class="out" markdown="1">
~~~
[moons 62e7791] Thinking about the moons
 2 files changed, 2 insertions(+)
 create mode 100644 moons.txt
~~~
</div>

Our repository is now in the state shown below:

<img src="img/git-branching-04.svg" alt="After Committing on Moons Branch" />

The `moons` branch has advanced to record the changes we just made,
but `master` is still where it was.
If we switch back to `master`:

<div class="in" markdown="1">
~~~
$ git checkout master
~~~
</div>

our changes seem to disappear:

<div class="in" markdown="1">
~~~
$ ls
~~~
</div>
<div class="out" markdown="1">
~~~
mars.txt
~~~
</div>
<div class="in" markdown="1">
~~~
$ cat mars.txt
~~~
</div>
<div class="out" markdown="1">
~~~
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~
</div>

They're still in the repository&mdash;they're just not in
the revision that `master` is currently pointing to.
In essence,
we've created a parallel timeline that shares some history with the original one before diverging.

Let's make some changes in the `master` branch to further illustrate this point:

<div class="in" markdown="1">
~~~
$ echo "Should we go with a classical name like Ares Base?" > names.txt
$ git status
~~~
</div>
<div class="out" markdown="1">
~~~
# On branch master
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#    names.txt
nothing added to commit but untracked files present (use "git add" to track)
~~~
</div>
<div class="in" markdown="1">
~~~
$ git add names.txt
$ git commit -m "We will need a cool name for our secret base"
~~~
</div>
<div class="out" markdown="1">
~~~
[master dfcf908] We will need a cool name for our secret base
 1 file changed, 1 insertion(+)
 create mode 100644 names.txt
~~~
</div>

Our repository is now in this state:

<img src="img/git-branching-05.svg" alt="After Committing on Master Branch" />

`master` and `moons` have both moved on from their original common state,
but in different ways.
They could continue independent existence indefinitely,
but at some point we'll probably want to [merge](../../gloss.html#merge) our changes.
Let's do that now:

<div class="in" markdown="1">
~~~
$ git branch
~~~
</div>
<div class="out" markdown="1">
~~~
* master
  moons
~~~
</div>
<div class="in" markdown="1">
~~~
$ git merge moons
~~~
</div>

When we run the `git merge` command,
Git opens an editor to let us write a log entry about what we're doing.
The editor session initially contains this:

<div class="file" markdown="1">
~~~
Merge branch 'moons'

# Please enter a commit message to explain why this merge is necessary,
# especially if it merges an updated upstream into a topic branch.
#
# Lines starting with '#' will be ignored, and an empty message aborts
# the commit.
~~~
</div>

If we notice that something is wrong and decide not to complete the merge,
we must delete everything in the file&mdash;Git interprets an empty log message to mean,
"Don't proceed."
Otherwise,
everything that isn't marked as a comment with `#` will be saved to the log.
In this case,
we'll stick with the default log message.
When we save the file and exit the editor, Git displays this:

<div class="out" markdown="1">
~~~
Merge made by the 'recursive' strategy.
 mars.txt  | 1 +
 moons.txt | 1 +
 2 files changed, 2 insertions(+)
 create mode 100644 moons.txt
~~~
</div>

We now have all of our changes in one place:

<div class="in" markdown="1">
~~~
$ ls
~~~
</div>
<div class="out" markdown="1">
~~~
mars.txt    moons.txt    names.txt
~~~
</div>

and our repository looks like this:

<img src="img/git-branching-06.svg" alt="After Merging" />

We can ask Git to draw a diagram of the repository's history with this command:

<div class="in" markdown="1">
~~~
$ git log --oneline --topo-order --graph
~~~
</div>
<div class="out" markdown="1">
~~~
*   e0cf8ab Merge branch 'moons'
|\
| * 62e7791 Thinking about the moons
* | dfcf908 We will need a cool name for our secret base
|/
* 005937f Thoughts about the climate
* 34961b1 Concerns about Mars's moons on my furry friend
* f22b25e Starting to think about Mars
~~~
</div>

This ASCII art is fine for small sets of changes,
but for anything significant,
it's much better to use a GUI
that can draw graphs using lines instead of characters.

Branching strikes most newcomers as unnecessary complexity,
particularly for single-author projects.
After all,
if we need to make some changes to a project,
what do we gain by creating parallel universes?

The answer is that branching makes it easy for us to concentrate on one thing at a time.
Suppose we are part-way through rewriting a function that calculates spatial correlations
when we realize that the task would be easier if our file I/O routines always stored things as complex numbers.
Most people would put the spatial correlation changes aside,
change the file I/O,
then (hopefully) come back to the original task.

The problem with this is that we have to remember what we were doing,
even if we realize halfway through rewriting file I/O that we should also rewrite our error handling.
It's quite common to wind up with half a dozen tasks stacked on top of one another,
and quite hard to them all straight.
Branching allows us to put what we're doing in a safe place,
solve the new problem,
then resume our earlier work.

In practice,
most developers never make changes directly in the `master` branch.
Instead,
they create a new branch from it for every change they want to make,
then merge those branches back to `master` when the work is complete.
Returning to our hypothetical example,
we would:

1.  create a branch called something like `better-spatial-correlation` for those changes;
2.  go back to master and create another branch called `file-input-produces-complex-values` for *those* changes;
3.  merge `file-input-produces-complex-values` into `master`;
4.  merge `master` into `better-spatial-correlation`; and
5.  finish work on the spatial correlation function and merge it all back into `master`.

And if,
partway through this process,
our supervisor asked us to change the graph-plotting routines to conform to the university's new style guide,
we would simply switch back to `master`,
create a branch for that,
make the changes,
produce the desired graphs,
and leave the changes parked in that branch until we were ready to merge them.
