---
layout: lesson
root: ../..
title: The Git Development Workflow
level: intermediate
---

An important best practice for scientific software developers is to
focus on making *incremental changes*.  In this section, we'll focus
on the basic techniques for saving your work using Git.  In this
lesson, you'll learn a third use for the checkout command, how to work
with development branches, and then how to save revisions in the
staging area.

## Objectives

After this lesson, students should be able to:

* Create and manipulate development branches using `git checkout` and
  `git branch`.
* Inspect the changes they have made to their files using `git status`
  and `git diff`.
* Add and remove changes from the staging area using `git add`, `git
  reset`, and `git rm`.
* Save versions of their work using `git commit`.

## Development Branches

In the last lesson we cloned a repository and explored its history.
We also learned about HEAD, a special *pointer* to the revision of our
project that we most recently checked out.  We call lightweight
pointers to revisions in Git *branches*.  This is confusing at
first, because many newcomers to Git think of branches as part of the
structure of the repository's history, instead of pointers to it.

The best way to think of branches is to first imagine your repository as a
notebook shared by several collaborators in your laboratory, with
bookmarks (the *branches*) on the page of their most recent experiment.
Every time somebody performs an experiment, they move their bookmark
to the page of the newest experiment.  Since laboratory paper is
cheap, researchers take notes of experiments that might fail, and
*reset* their notes whenever they don't like an outcome.  Any time
somebody wants to keep track of one of their collaborators, they just
need to look up their bookmarks.  Since bookmarks are cheap, some of
the scientists track multiple experiments, even risky ones that might
not pan out, with their bookmarks.  The bookmarks themselves are not
the entries in the notebook, but they are a handy way to keep track of
active work.

In Git, a branch is a reference that follows your development.  When
you create new revisions, your current branch will automatically move
to new commits you make.  Branches are also primary tools for sharing
and collaboration, as they allow you to refer others to
work-in-progress (and scientific software development is always
work-in-progress).

When a repository is created, Git creates a branch called `master`.
Unlike HEAD, the only thing special about `master` is that it's the
default name given to this branch.  It's generally not a great idea to
work directly on `master`, so we'll start this lesson by creating a new
branch from it for development.

By default, Git creates branches that point to the currently
checked-out revision.  We would like to start our development at the
most up-to-date version of the repository, so we check out our local
copy of the `master` branch first.  This `checkout` command performs
two actions.  First, it checks out the repository to the revision
pointed to by the `master` branch, in this case: `61fd2bc`.  Second,
it *activates* the `master` branch.  Keep in mind, Git will
automatically move the active branch pointer for you when you make
commits.

~~~
$ git checkout master
Switched to branch 'master'
~~~

We can create and activate a branch in the same step by adding the
`-b` flag to `git checkout`.  Since later we'll be fixing an error
in the  file `python_pipeline.ipy`, we'll name the branch
`fix_pipeline`.

~~~
$ git checkout -b pipeline_fix
Switched to a new branch 'pipeline_fix'
~~~

We can see the current branches available with `git branch`, but this
command is much more useful if we add the `-v` flag (for verbose).
Note that `git branch` lists the local branches in our repository,
`master` and `pipeline_fix`, which we just created.  The `*` next to
`pipeline_fix` indicates that it is active.

~~~
$ git branch -v
  master       61fd2bc Made fixes to Python pipeline
* pipeline_fix 61fd2bc Made fixes to Python pipeline
~~~

Now that we've created our development branch, let's start preparing
some changes to our repository.

### Checkpoint 1

* **[1A]** Why do `master` and `pipeline_fix` point to the same
  revision even though they are different branches?
* **[1B]** The `--decorate` argument to `git log` will annotate your
  history view with the location of your current branches.  Try
  calling `git log --oneline --decorate`.  Discuss with your partner
  what you think `origin/master` refers to.
* **[1C]** The `DETACHED HEAD` state occurs when the repository does
  not have an active branch, frequently after checking out a commit or
  remote reference.  Try calling the command `git checkout 61fd` and
  explaining the output message to your partner.
* **[1D]** You can use `git branch -D` to delete a branch that is not
  active.  Remember that this only deletes your local reference to a
  commit.  Try deleting the `master` branch in your repository by
  calling `git branch -D master`.


## Seeing Changes

In the previous lesson, we used `git checkout` to undo some changes
we had made to a file.  We also learned about `diff` output, which is
a common way for computers to tell us how two text files differ.  In
this section, we'll learn about the `git status` and `git diff`
commands, which help us see changes in our repository.  We'll also learn
how to use the `git show` command to retrieve a copy of a file from
anywhere in our repository's history.

Let's start by repeating our dangerous idea from the previous lesson,
and delete an important data file from the repository!

~~~
$ rm Lumi.2763.csv
$ ls Lumi.2763.csv                                                                            ✖
ls: Lumi.2763.csv: No such file or directory
~~~

Now that the file is gone, we'd like to ask Git if it was paying
attention and knows something is missing.  We're going to do that now
with one of the most powerful tools in the Git toolbox, the `git
status` command.

~~~
$ git status
On branch pipeline_fix
Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    Lumi.2763.csv

no changes added to commit (use "git add" and/or "git commit -a")
~~~

Whew, that's a lot of output.  Let's take it line by line.  First,
the `git status` command reports that we're on the `pipeline_fix`
branch (if it doesn't, and you don't understand why, re-read the
previous section).  Second, Git is reporting that you have made some
changes, but these are not yet saved.  Git refers to your files as the
*working directory*, or sometimes as the *working tree*.  If you see
*working* in the name, be confident that Git is referring to the files
you can see in front of you with commands like `ls`.  Finally, Git
summarizes the changes.  In this case, since a file was completely
deleted, it reports that `Lumi.2763.csv` has been deleted.  Let's
ignore the deleted file for now, and make some changes to
`python_pipeline.ipy`.  I encourage you to make real changes with your
editor, but for the purposes of this example, I'm going to append a
line of text to the end of `python_pipeline.ipy` from the command line.

~~~
$ echo "# This is a comment." >> python_pipeline.ipy
~~~

We call `git status` to see the effects of modifying
`python_pipeline.ipy`.

~~~
$ git status
On branch pipeline_fix
Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    Lumi.2763.csv
	modified:   python_pipeline.ipy

no changes added to commit (use "git add" and/or "git commit -a")
~~~

Okay!  Now we're going to show you one other type of change you might
make to a repository, creating a new file.  Let's call this file,
`documentation.txt`, and let's put the line: `TODO` in by itself.
Again, please use your favorite editor for this, here's a command line
script for the purposes of this example.

~~~
echo "TODO" > documentation.txt
~~~

We're going to call `git status` one more time.

~~~
$ git status
On branch pipeline_fix
Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    Lumi.2763.csv
	modified:   python_pipeline.ipy

Untracked files:
  (use "git add <file>..." to include in what will be committed)

	documentation.txt

no changes added to commit (use "git add" and/or "git commit -a")
~~~

When we deleted and modified files already in the repository, `git
status` reported these as *unstaged* changes.  New files, however, are
*untracked*.  The difference is subtle.  New files have no Git
history, so Git ignores them by default (and tells you it is doing
so).

Sometimes it isn't enough to see that files have changed.  We want to
know how they changed!  The `git diff` command provides a report of
the *differences* in the working directory:

~~~
$ git diff
diff --git a/Lumi.2763.csv b/Lumi.2763.csv
deleted file mode 100755
index fb57a5b..0000000
--- a/Lumi.2763.csv
+++ /dev/null
@@ -1,10 +0,0 @@
-    ,  1   ,  2   ,  3   ,  4   ,  5   ,  6   ,  7   ,  8   ,  9   ,  10  ,  11  ,  12
-  A , 8053 , 9397 ,13857 ,  886 ,  985 , 1204 ,34036 ,34673 ,40062 , 1354 ,  944 , 1014
-  B , 6538 , 7707 , 7419 , 1098 , 1222 , 1221 ,30437 ,32832 ,31382 , 1217 ,  987 ,  967
-  C , 4734 , 4971 , 4376 , 1205 , 1386 , 1369 ,26260 ,24249 ,21131 , 1101 ,  871 ,  826
-  D , 3785 , 2310 , 3154 , 1088 , 1085 , 1038 ,14243 ,12519 ,13091 ,  699 ,  822 ,  772
-  E , 2538 , 1728 , 2338 , 1063 , 1182 , 1202 ,10260 , 6621 , 6756 ,  922 ,  857 ,  731
-  F , 1453 , 1128 , 1255 , 1022 ,  951 , 1302 , 2404 , 2900 , 2161 ,  760 ,  640 ,  617
-  G , 1409 , 1010 , 1080 , 1019 , 1309 ,  975 , 1943 , 1627 , 1495 ,  737 ,  701 ,  696
-  H , 1476 , 1025 , 1164 , 1120 , 1312 , 1310 , 1916 , 1559 , 1711 , 1040 ,  888 ,  927
-
diff --git a/python_pipeline.ipy b/python_pipeline.ipy
index 2ecc443..1c43ea6 100644
--- a/python_pipeline.ipy
+++ b/python_pipeline.ipy
@@ -43,3 +43,4 @@ colorbar()
 norm_data = h/(bgal_data-0.07)
 imshow(norm_data)

+# This is a comment.
~~~

We discussed `diff` output in the previous lesson.  Here `git diff`
is showing us the changes in our
current working directory.  Again, this is a lot of output, but for
the most part, we only pay attention to the lines that begin with a
`-` or a `+`.  We see that all of the lines in `Lumi.2763.csv` have
been deleted, which is what we would expect for a deleted file.  In
`python_pipeline.ipy`, our line `#This is a comment.` has been added.
The content of our new file is not shown.  Why is that?  Git assumes
that you only want to see the changes to files you've explicitly told
it about.  How do you tell Git that you'd like to track a new file?
We do this with the `git add` command, which we'll cover in the next
section.

### Checkpoint 2

* **[2A]** What is the difference between the way Git uses the words
  *unstaged* and *untracked*?  Try explaining this to your partner.
* **[2B]** The `git status` and `git diff` commands take a few useful
  flags.  Try calling `git status -s` and `git diff --stat` to see
  terse summaries of the status of the repository and changes that
  have been made.
* **[2C]** Sometimes it is useful to see the state of a file before you
  started modifying it.  The `git show` command provides this
  capability, although the syntax is slightly awkward.  Try calling
  `git show HEAD:Lumi.2763.csv` and `git show master:Lumi.2763.csv`,
  then explain the difference between the two *revision parameters*,
  `master`, and `HEAD`.  Consult
  [the man page for git revisions](https://www.kernel.org/pub/software/scm/git/docs/gitrevisions.html)
  for some gory details.


## Staging Changes

So far, we have made three changes to the working directory.  We
deleted a file, we modified a file, and we created a new one.  Let's
suppose that we want to leave the deleted file alone for now, and save
the other changes to our history.  Since we frequently make changes
across many of our files while we are working, but are only interested
in saving some of these changes into each revision of the repository,
Git introduces the concept of *staging*.

To communicate which changes we would like to save, we put them in the
*stage*.  Think of the stage like a set at a camera studio.  We can
add and remove files from the stage in preparation for saving our
changes as a revision in history. Everything on the stage will be in
the picture when it's taken.  Keep in mind that Git-speak for "create
a new revision with these changes I have staged" is *commit*.
The word "commit" is shamelessly used both as a verb for referring to the
process of creating revisions, and as a noun for referring to a
specific revision.

*Hint: You may see the words index and cache in the documentation or
 while reading about Git.  These terms both refer to the stage!*

Let's get our keyboards dirty!  To add files to the stage, we use the
`git add` command.  Let's try that now with our
untracked-until-right-now `documentation.txt`.

~~~
$ git add documentation.txt
~~~

Let's see what's changed in `git status`:

~~~
$ git status
On branch pipeline_fix
Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

	new file:   documentation.txt

Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    Lumi.2763.csv
	modified:   python_pipeline.ipy
~~~

`documentation.txt` has graduated from an untracked file to a change
ready to be saved in our next revision!  We can also use `git add` to
stage changes to modified files.

~~~
$ git add python_pipeline.ipy
$ git status
On branch pipeline_fix
Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

	new file:   documentation.txt
	modified:   python_pipeline.ipy

Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    Lumi.2763.csv
~~~

Let's see how this affected the output of `git diff`.

~~~
$ git diff
diff --git a/Lumi.2763.csv b/Lumi.2763.csv
deleted file mode 100755
index fb57a5b..0000000
--- a/Lumi.2763.csv
+++ /dev/null
@@ -1,10 +0,0 @@
-    ,  1   ,  2   ,  3   ,  4   ,  5   ,  6   ,  7   ,  8   ,  9   ,  10  ,  11  ,  12
-  A , 8053 , 9397 ,13857 ,  886 ,  985 , 1204 ,34036 ,34673 ,40062 , 1354 ,  944 , 1014
-  B , 6538 , 7707 , 7419 , 1098 , 1222 , 1221 ,30437 ,32832 ,31382 , 1217 ,  987 ,  967
-  C , 4734 , 4971 , 4376 , 1205 , 1386 , 1369 ,26260 ,24249 ,21131 , 1101 ,  871 ,  826
-  D , 3785 , 2310 , 3154 , 1088 , 1085 , 1038 ,14243 ,12519 ,13091 ,  699 ,  822 ,  772
-  E , 2538 , 1728 , 2338 , 1063 , 1182 , 1202 ,10260 , 6621 , 6756 ,  922 ,  857 ,  731
-  F , 1453 , 1128 , 1255 , 1022 ,  951 , 1302 , 2404 , 2900 , 2161 ,  760 ,  640 ,  617
-  G , 1409 , 1010 , 1080 , 1019 , 1309 ,  975 , 1943 , 1627 , 1495 ,  737 ,  701 ,  696
-  H , 1476 , 1025 , 1164 , 1120 , 1312 , 1310 , 1916 , 1559 , 1711 , 1040 ,  888 ,  927
-
~~~

Why do our staged changes no longer appear in the output of `git
diff`?  The answer to this question is subtle.  When we change
revisions using `git checkout`, Git sets the contents of our working
directory *and* the staging area to the content of the revision.  When
we call `git diff`, it is comparing our working directory to the
staging area.  If we want to see how our *staged* changes will look in
history after they've been commited, we add the `--staged` flag to the
`git diff` command.

~~~
$ git diff --staged
diff --git a/documentation.txt b/documentation.txt
new file mode 100644
index 0000000..1333ed7
--- /dev/null
+++ b/documentation.txt
@@ -0,0 +1 @@
+TODO
diff --git a/python_pipeline.ipy b/python_pipeline.ipy
index 2ecc443..1c43ea6 100644
--- a/python_pipeline.ipy
+++ b/python_pipeline.ipy
@@ -43,3 +43,4 @@ colorbar()
 norm_data = h/(bgal_data-0.07)
 imshow(norm_data)

+# This is a comment.
~~~

Sometimes we make a mistake and stage some changes we didn't mean to.
It would be easy to remember how to undo changes if Git just used the
commands `stage` and `unstage`.  Unfortunately, Linus Torvalds gave us
the confusing pair: `add` and `reset`.  If you think of the
photography studio metaphor, then *resetting* the stage to what it
looked like in the last photo kind of works, but it's still confusing.

Let's reset the staged version of `python_pipeline.ipy`.  This command
only affects the staging area, so the file will remain changed, but
the changes will go from staged to unstaged:

~~~
$ git reset python_pipeline.ipy
Unstaged changes after reset:
D	Lumi.2763.csv
M	python_pipeline.ipy
$ git status
On branch pipeline_fix
Changes to be committed:
  (use "git reset HEAD <file>..." to unstage)

	new file:   documentation.txt

Changes not staged for commit:
  (use "git add/rm <file>..." to update what will be committed)
  (use "git checkout -- <file>..." to discard changes in working directory)

	deleted:    Lumi.2763.csv
	modified:   python_pipeline.ipy
~~~


### Checkpoint 3

* **[3A]** The `--` flag is used to separate filenames and paths from
  the other parameters in a `git` command line argument.  What is the
  output of `git diff -- python_pipeline.ipy`?  What other commands
  take a single filename or set of files as an argument?
* **[3B]** Git doesn't let you stage the removal of a file with the
  `git add` command.  Try using the command `git rm Lumi.2763.csv`.
  To unstage this removal, the appropriate command is `git reset --
  Lumi.2763.csv`.  How would you restore the last revision of
  `Lumi.2763.csv` from the repository's history?
* **[3C]** Use your editor to make changes to files in the
  repository.  Experiment with staging changes with `git add`,
  unstaging them with `git reset`, deleting files with `git rm`,
  restoring them with `git checkout`, and visualizing the staging area
  with `git status` and `git diff`.


## Committing Changes

* Walk through `git commit` using the `nano` editor.
* Use email metaphor for `git commit`.
* make a couple of commits on the development branch.
* use `git log --decorate` to see revisions, how they differ.
* Explain how Git commits are like lightweight publishing.
  + Simple to make, yet allow you to easily "sign" revisions of very
	complicated objects.
  + Like a lab notebook.  Can be well-categorized or it can be very
	sketchy.  Use it the way that it complements your workflow.  But
	always commit at the end of the day, even if it's something you
	will throw away later.

### Integrative Exercise

* Checkout master.  Make a commit.  Now swap the positions of the two
  branches, by either moving the branches with `git reset` or by
  deleting them and recreating them in the appropriate positions.

## Review

* branches
* viewing changes
* staging changes
* committing changes
