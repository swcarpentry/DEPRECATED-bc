---
layout: lesson
root: ../..
title: A Better Kind of Backup
level: novice
---
The first time we use Git on a new machine,
we need to configure a few things (we'll insert blank lines
between groups of shell commands to make them easier to read):

~~~
$ git config --global user.name "Vlad Dracula"

$ git config --global user.email "vlad@tran.sylvan.ia"

$ git config --global color.ui "auto"

$ git config --global core.editor "nano"
~~~

(Please use your own name and email address instead of Dracula's,
and please make sure you choose an editor that's actually on your system
rather than `nano`.)

Git commands are written `git verb`,
where `verb` is what we actually want it to do.
In this case,
we're telling Git:

*   our name and email address,
*   to colorize output,
*   what our favorite text editor is, and
*   that we want to use these settings globally (i.e., for every project),

The four commands above only need to be run once:
Git will remember the settings until we change them.
Once Git is configured,
we can start using Git.
Let's create a directory for our work:

~~~
$ mkdir planets

$ cd planets
~~~

and tell Git to make it a [repository](../gloss.html#repository):

~~~
$ git init
~~~

If we use `ls` to show the directory's contents,
it appears that nothing has changed:

~~~
$ ls
~~~

But if we add the `-a` flag to show everything,
we can see that Git has created a hidden directory called `.git`:

~~~
$ ls -a
.	..	.git
~~~

Git stores information about the project in this special sub-directory.
If we ever delete it,
we will lose the project's history.

We can ask Git for the status of our project at any time like this:

~~~
$ git status
# On branch master
#
# Initial commit
#
nothing to commit (create/copy files and use "git add" to track)
~~~

We'll explain what `branch master` means later.
For the moment,
let's add some notes about Mars's suitability as a base.
(We'll `cat` the text in the file after we edit it so that you can see what we're doing,
but in real life this isn't necessary.)

~~~
$ nano mars.txt

$ cat mars.txt
Cold and dry, but everything is my favorite color

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
~~~

In order to understand what the "untracked files" message really means, we need to understand the basic model of how git groups your files.  There are essentially three places your work can reside:
- As a 'revision'.  A revision is like a snapshot of your work which you can come back to at any time.  Revisions also have a relational tree-like structure; we'll get into the details of this later, but for now just imagine a very simple chain of revisions; every time you make a new revision, it has (for now) one parent revision, which is just whatever revision came before it.
- In the 'index'.  Git's index is like the staging area where you put your work to tell git 'these are the changes I'd like to package as my next revision'.
- In the 'working tree'.  Git's working tree is just the new work you've done that hasn't been added to the index or committed; every time you save a change to a file in the way you're used to, you've put new stuff in the working tree.  

With that mental model in mind, the "untracked files" message means that there's a file in the directory
that Git isn't keeping track of - it's in your working tree, but not your index.
We can tell Git that it should add it to the index so that we can later commit it:

~~~
$ git add mars.txt
~~~

and check that the right thing happened like this:

~~~
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
~~~

mars.txt is now in the index - Git now knows that it's supposed to keep track of this file,
but it hasn't yet recorded any changes for posterity as a commit.
To get it to do that,
we need to run one more command:

~~~
$ git commit -m "Starting to think about Mars"
[master (root-commit) f22b25e] Starting to think about Mars
 1 file changed, 1 insertion(+)
 create mode 100644 mars.txt
~~~

When we run `git commit`,
Git takes everything we have told it to save by using `git add`
and stores a copy permanently inside the special `.git` directory.
This permanent copy is called a [revision](../../gloss.html#revision). 
In git parlance, 'git add mars.txt' promoted mars.txt from the working tree to the index,
and 'git commit -m ....' packaged up everything in the index as a commit that we can return to at a later time if need be. 
We use the `-m` flag (for "message")
to record a comment that will help us remember later on what we did and why.
If we just run `git commit` without the `-m` option,
Git will launch `nano` (or whatever other editor we configured at the start)
so that we can write a longer message.

If we run `git status` now:

~~~
$ git status
# On branch master
nothing to commit, working directory clean
~~~

it tells us everything is up to date.
If we want to know what we've done recently,
we can ask Git to show us the project's history:

~~~
$ git log
commit f22b25e3233b4645dabd0d81e651fe074bd8e73b
Author: Vlad Dracula <vlad@tran.sylvan.ia>
Date:   Thu Aug 22 09:51:46 2013 -0400

    Starting to think about Mars
~~~

Now suppose Dracula adds more information to the file:

~~~
$ nano mars.txt

$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
~~~

If we run `git status`,
it tells us the file has been modified:

~~~
$ git status
# On branch master
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
no changes added to commit (use "git add" and/or "git commit -a")
~~~

The last line is the key phrase:
"no changes added to commit".
We have changed this file in our working tree, but we haven't promoted those changes to the index or saved them as as commit. 
Let's double-check our work using `git diff`,
which shows us the differences between
the current state of the file
and the most recently saved version:

~~~
$ git diff
diff --git a/mars.txt b/mars.txt
index df0654a..315bf3a 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,2 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
~~~

The output is cryptic because
it is actually a series of commands for tools like editors and `patch`
telling them how to reconstruct one file given the other.
If we can break it down into pieces:

1.   The first line tells us that Git is using the Unix `diff` command
     to compare the old and new versions of the file.
2.   The second line tells exactly which [revisions](../../gloss.html#revision) of the file
     Git is comparing;
     `df0654a` and `315bf3a` are unique computer-generated labels for those revisions.
3.   The remaining lines show us the actual differences
     and the lines on which they occur.
     The numbers between the `@@` markers indicate which lines we're changing;
     the `+` on the lines below show that we are adding lines.

Let's commit our change:

~~~
$ git commit -m "Concerns about Mars's moons on my furry friend"
# On branch master
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
no changes added to commit (use "git add" and/or "git commit -a")
~~~

Whoops:
Git won't commit because we didn't use `git add` first - there's nothing in the index and nothing for git to make a commit out of!
Remember to promote our work from the working tree to the index first using 'git add':

~~~
$ git add mars.txt

$ git commit -m "Concerns about Mars's moons on my furry friend"
[master 34961b1] Concerns about Mars's moons on my furry friend
 1 file changed, 1 insertion(+)
~~~

Git insists that we add files to the set we want to commit
before actually committing anything
because we often won't commit everything at once.
For example,
suppose we're adding a few citations to our supervisor's work
to our thesis.
We might want to commit those additions,
and the corresponding addition to the bibliography,
but *not* commit the work we've been doing on the conclusion.
To allow for this,
Git has a special staging area
where it keeps track of things that have been added to
the current [change set](../gloss.html#change-set)
but not yet committed.
`git add` puts things in this area (the index),
and `git commit` then copies them to long-term storage (as a commit):

<img src="img/git-staging-area.svg" alt="The Git Staging Area" />

The following commands show this in action:

~~~
$ nano mars.txt

$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity

$ git diff
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~

So far, so good:
we've made a change,
and `git diff` tells us what it is.
Now let's put that change in the staging area
and see what `git diff` reports:

~~~
$ git add mars.txt

$ git diff
~~~

There is no output:
as far as Git can tell,
there's no difference between what it's been asked to save permanently
and what's currently in the directory.
However,
if we do this:

~~~
$ git diff --staged
diff --git a/mars.txt b/mars.txt
index 315bf3a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~

it shows us the difference between
the last committed change
and what's in the staging area.
Let's save our changes:

~~~
$ git commit -m "Thoughts about the climate"
[master 005937f] Thoughts about the climate
 1 file changed, 1 insertion(+)
~~~

check our status:

~~~
$ git status
# On branch master
nothing to commit, working directory clean
~~~

and look at the history of what we've done so far:

~~~
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
~~~

If we want to see what we changed when,
we use `git diff` again,
but refer to old versions
using the notation `HEAD~1`, `HEAD~2`, and so on:

~~~
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
~~~

Recall above we mentioned that revisions have a relational structure, for now just like a simple chain; in git, the word `HEAD` always refers to the most recent end of that chain, the last revision you tacked on.  In other words, `HEAD` means "the most recently saved version".  Every time you do git commit, a new revision is tacked onto the end of that chain, and `HEAD` moves forward to point at that new latest revision.  We can step backwards on the chain using the `~` notation;
`HEAD~1` (pronounced "head minus one")
means "the previous revision", or `HEAD~12` means "12 revisions ago".
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

~~~
$ git diff f22b25e3233b4645dabd0d81e651fe074bd8e73b mars.txt
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~

That's the right answer,
but typing random 40-character strings is annoying,
so Git lets us use just the first few:

~~~
$ git diff f22b25e mars.txt
diff --git a/mars.txt b/mars.txt
index df0654a..b36abfd 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~

All right:
we can save changes to files and see what we've changed---how
can we restore older versions of things?
Let's suppose we accidentally overwrite our file:

~~~
$ nano mars.txt

$ cat mars.txt
We will need to manufacture our own oxygen
~~~

`git status` now tells us that the file has been changed,
but those changes haven't been staged:

~~~
$ git status
# On branch master
# Changes not staged for commit:
#   (use "git add <file>..." to update what will be committed)
#   (use "git checkout -- <file>..." to discard changes in working directory)
#
#	modified:   mars.txt
#
no changes added to commit (use "git add" and/or "git commit -a")
~~~

We can put things back the way they were like this:

~~~
$ git reset --hard HEAD
HEAD is now at 005937f Thoughts about the climate

$ cat mars.txt
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~

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

In other words, `git reset --hard HEAD` wipes out our index and makes our working tree match `HEAD` exactly - so watch out!  Any work you've done since your last commit is now gone forever.

There are less drastic things we can do with git reset as well:  

```
$ git reset --mixed HEAD
```

will empty the index just like `git reset --hard`, but it won't touch your working tree, so you won't lose any work (much safer) - use this if you've added some things with `git add` which you regret, and would like to empty the index (but not the working tree).  Even more tame is the ever-gentle

```
$ git reset --soft HEAD~
```

which will leave your working tree alone like `--mixed`, but also automatically stage all the differences between `HEAD` and `HEAD~` in the index, ready to be committed as a descendent of `HEAD~` - use this to change the commit message or add something you forgot to the last commit.

We can summarize the impact of the three options as follows...

| Command             | Effect                                         | Common Use                                 |
|---------------------|------------------------------------------------|--------------------------------------------|
| `git reset --soft`  | No changes to working tree or index            |Change the commit message you just made     |
| `git reset --mixed` | Remove staged changes from index               |Change what's going into the next commit    |
| `git reset --hard`  | Remove changed files in index and working tree |Completely abandon everything and go back   |


In all cases, `git reset` is moving where `HEAD` is currently pointing; think of `HEAD` as the end of that chain of commits that git is making, and moving it changes where the next commit is going to get tacked on to; the difference between `--hard`, `--mixed` and `--soft` is just what happens to the index and working trees as we move `HEAD`.  In later lessons we'll learn more about the tremendous power of git's relational structure, but for now it's enough to understand that making a new commit moves the `HEAD` forward by one commit, and reset moves `HEAD` backwards along the chain to wherever we reset to.


But what if we want to recover somes files without wandering around on our chain of commits and messing with where `HEAD` is?
For example,
what if we have added some material to the conclusion of our paper that we'd like to keep,
but we want to get back an earlier version of the introduction, all without touching our sequence of revisions?
In that case,
we want to check out an older revision of the file,
so we do something like this:

~~~
$ git checkout 123456 mars.txt
~~~

but use the first few digits of an actual revision number instead of 123456.
To get the right answer,
we must use the revision number that identifies the state of the repository
*before* the change we're trying to undo.
A common mistake is to use the revision number of
the commit in which we made the change we're trying to get rid of:

<img src="img/git-when-revisions-updated.svg" alt="When Git Updates Revision Numbers" />

The fact that files can be reverted one by one
tends to change the way people organize their work.
If everything is in one large document,
it's hard (but not impossible) to undo changes to the introduction
without also undoing changes made later to the conclusion.
If the introduction and conclusion are stored in separate files,
on the other hand,
moving backward and forward in time becomes much easier.


What if we have files that we do not want Git to track for us? For
example, we might have backup files created by the editor, or result
files or intermediate data files that we do not need to track, and
these can clutter up our `git status` output.

Let's create some data files that we don't want to track:

```
$ touch a.dat b.dat c.dat
$ mkdir results
$ touch results/a.out results/b.out
```

Then, we can look at the status:

```
$ git status
# On branch master 
# 
# Initial commit 
# 
# Untracked files: 
# (use "git add <file>..." to include in what will be committed) 
# 
# a.dat 
# b.dat 
# c.dat 
# results/ 
nothing added to commit but untracked files present (use "git add" to track)
```

Seeing all of these files can be annoying, and worse, the clutter can
hide important files from us. Let's tell Git to ignore those files.

Git has a powerful mechanism for ignoring files.  We will look at a
basic example.  In a repository, the file `.gitignore` is used to
specify files to ignore.

```
$ nano .gitignore
$ cat .gitignore
*.dat
/results/
```

This tells Git to ignore any file ending with `.dat`, and all files
under the `results` directory. (If any of these files were already being
tracked, they would not be affected by `.gitginore`.)

Now, we have a much cleaner status.


```
# On branch master
#
# Initial commit
#
# Untracked files:
#   (use "git add <file>..." to include in what will be committed)
#
#	.gitignore
nothing added to commit but untracked files present (use "git add" to track)
```


Now Git notices the `.gitignore` file.  You might think we wouldn't
want to track this file, but if our programs are creating these extra
data files, all users of this repository will probably want to ignore
them too.  So, we can add the ignore file to the repository and check
it in.

```
$ git add .gitignore
$ git commit -m "Add the ignore file"
$ git status
# On branch master
nothing to commit, working directory clean
```

When `git status` is clean, it is much easier to tell that all of our
work is checked in.

We can see the ignored files in the status if we want.  That can be
helpful when making sure that no files were accidentally ignored.

```
$ git status --ignored
# On branch master
# Ignored files:
#  (use "git add -f <file>..." to include in what will be committed)
#
#        a.dat
#        b.dat
#        c.dat
#        results/

nothing to commit, working directory clean
```

Finally, the ignore file helps us avoid accidentally adding files to
the repository that we don't want.

```
$ git add a.dat
The following paths are ignored by one of your .gitignore files:
a.dat
Use -f if you really want to add them.
fatal: no files added
```

Git notices that `a.dat` would be excluded by the rules given in the
`.gitignore` file and doesn't automatically include it. We can
override this by using the `-f` flag.

As shown here, the `.gitignore` file just works on a specific
repository, but you can also set global rules for `git ignore` using a
`.gitginore` file in your home directory.
