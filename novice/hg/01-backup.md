---
layout: lesson
root: ../..
title: A Better Kind of Backup
---
<div class="objectives" markdown="1">

#### Objectives
* Explain which initialization and configuration steps are required once per machine,
  and which are required once per repository.
* Add files to Mercurial's collection of tracked files.
* Go through the modify-commit cycle for single and multiple files
  and explain where information is stored before and after the commit.
* Identify and use Mercurial revision numbers and changeset identifiers.
* Compare files with old versions of themselves.
* Restore old versions of files.
* Configure Mercurial to ignore specific files,
  and explain why it is sometimes useful to do so.

</div>

We'll start by exploring how version control can be used
to keep track of what one person did and when.
Even if you aren't collaborating with other people,
version control is much better for that than this:

<div>
  <a href="http://www.phdcomics.com"><img src="img/phd101212s.gif" alt="Piled Higher and Deeper by Jorge Cham, http://www.phdcomics.com" /></a>
  <p>"Piled Higher and Deeper" by Jorge Cham, http://www.phdcomics.com</p>
</div>

#### Setting Up

The first time we use Mercurial on a new machine,
we need to configure a few things.

Dracula sets up his new Windows laptop by using his editor to create
a new file called `%USERPROFILE%\Mercurial.ini`
(that's spelled `$USERPROFILE/Mercurial.ini` if you are in `gitbash`)
containing the following
lines:

~~~
[ui]
username = Vlad Dracula <vlad@tran.sylvan.ia>
editor = nano

[extensions]
color =

[color]
mode = win32
~~~

Wolfman has both a Mac laptop and a Linux one and he uses his editor
to create a new file called `~/.hgrc` on both of those machines with
the same contents:

~~~
[ui]
username = Jack Wolfman <jack@cali.forn.ia>
editor = nano

[extensions]
color =
~~~

(Please use your own name and email address instead of Dracula's
or Wolfman's,
and make sure you choose an editor that's actually on your system,
such as `notepad` on Windows.)

Those configuration file settings tell Mercurial:

* our name and email address,
* what our favorite text editor is, and
* to colorize output.

The fact that these settings are in the Mercurial configuration file in
our home directory means that they will be used for every project on this
machine.
This bit of setup only needs to be done once.

#### Creating a Repository

Once Mercurial is configured,
we can start using it.
Let's create a directory for our work:

~~~
$ mkdir planets
$ cd planets
~~~
{:class="in"}

and tell Mercurial to make it a [repository](../../gloss.html#repository)&mdash;a place where
Mercurial can store old versions of our files:

~~~
$ hg init
~~~
{:class="in"}

Mercurial commands are written `hg verb`,
where `verb` is what we actually want it to do.

If we use `ls` to show the directory's contents,
it appears that nothing has changed:

~~~
$ ls
~~~
{:class="in"}

But if we add the `-a` flag to show everything,
we can see that Mercurial has created a hidden directory called `.hg`:

~~~
$ ls -a
~~~
{:class="in"}
~~~
.	..	.hg
~~~
{:class="out"}

Mercurial stores information about the project in this special sub-directory.
If we ever delete it,
we will lose the project's history.

We can check that everything is set up correctly
by asking Mercurial to verify the structure of our repository:

~~~
$ hg verify
~~~
{:class="in"}
~~~
checking changesets
checking manifests
crosschecking files in changesets and manifests
checking files
0 files, 0 changesets, 0 total revisions
~~~
{:class="out"}

#### Tracking Changes to Files

Let's create a file called `mars.txt` that contains some notes
about the Red Planet's suitability as a base.
(We'll use `nano` to edit the file;
you can use whatever editor you like.
In particular, this does not have to be the editor that you set in your
Mercurial global configuration earlier.)

~~~
$ nano mars.txt
~~~
{:class="in"}

`mars.txt` has now been created and it contains a single line:

~~~
$ ls
~~~
{:class="in"}
~~~
mars.txt
~~~
{:class="out"}
~~~
$ cat mars.txt
~~~
{:class="in"}
~~~
Cold and dry, but everything is my favorite color
~~~
{:class="out"}


We can ask Mercurial to tell us what it knows about the files in our project
with the `hg status` command.
Mercurial tells us that it's noticed the new file:

~~~
$ hg status
~~~
{:class="in"}
~~~
? mars.txt
~~~
{:class="out"}

The `?` at the beginning of the line means that Mercurial isn't keeping
track of the file.
We can tell Mercurial that it should do so using `hg add`:

~~~
$ hg add mars.txt
~~~
{:class="in"}

and then check that the right thing happened:

~~~
$ hg status
~~~
{:class="in"}
~~~
A mars.txt
~~~
{:class="out"}

Mercurial now knows that it's supposed to keep track of `mars.txt`,
but it hasn't yet recorded any changes for posterity as a commit.
To get it to do that,
we need to run one more command:

~~~
$ hg commit -m "Starting to think about Mars"
~~~
{:class="in"}

When we run `hg commit`,
Mercurial takes the file we have told it about by using `hg add`
and stores a copy permanently inside the special `.hg` directory.

We use the `-m` flag (for "message")
to record a comment that will help us remember later on what we did and why.
If we just run `hg commit` without the `-m` option,
Mercurial will launch `nano` (or whatever other editor we configured at the start)
so that we can write a longer message.

If we run `hg status` now:

~~~
$ hg status
~~~
{:class="in"}

we get no output because everything is up to date.

If we want to know what we've done recently,
we can ask Mercurial to show us the project's history using `hg log`:

~~~
$ hg log
~~~
{:class="in"}
~~~
changeset:   0:72ab25fa99a1
tag:         tip
user:        Vlad Dracula <vlad@tran.sylvan.ia>
date:        Mon Apr 14 14:41:58 2014 -0400
summary:     Starting to think about Mars

~~~
{:class="out"}

`hg log` lists all changes committed to a repository
in reverse chronological order.
The listing for each [changeset](../../gloss.html#change-set) includes:

* the changeset's revision number and identifier
  (`0` and `72ab25fa99a1` in this case,
  but your identifier will likely be different),
* its tags
  (more about tags later),
* the changeset's author,
* when it was created,
* and the log message Mercurial was given when the changeset was created.

The revision number is a convenient integer shorthand for the hexidecimal
identifier.

> #### Where Are My Changes?
>
> If we run `ls` at this point, we will still see just one file called `mars.txt`.
> That's because Mercurial saves information about files' history
> in the special `.hg` directory mentioned earlier
> so that our filesystem doesn't become cluttered
> (and so that we can't accidentally edit or delete an old version).

#### Changing a File

Now suppose Dracula adds more information to the file.
(Again, we'll edit with `nano` and then `cat` the file to show its contents;
you may use a different editor, and don't need to `cat`.)

~~~
$ nano mars.txt
$ cat mars.txt
~~~
{:class="in"}
~~~
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
~~~
{:class="out"}

When we run `hg status` now,
it tells us that a file it already knows about has been modified:

~~~
$ hg status
~~~
{:class="in"}
~~~
M mars.txt
~~~
{:class="out"}

The `M` at the beginning of the line means that Mercurial has noticed that
we have modified the `mars.txt` file.

We can double-check our work using `hg diff`,
which shows us the differences between
the current state of the file
and the most recently committed version:

~~~
$ hg diff
~~~
{:class="in"}
~~~
diff -r 72ab25fa99a1 mars.txt
--- a/mars.txt  Mon Apr 14 14:41:58 2014 -0400
+++ b/mars.txt  Mon Apr 14 15:48:53 2014 -0400
@@ -1,1 +1,2 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
~~~
{:class="out"}

The output is cryptic because
it is actually a series of commands for tools like editors and `patch`
telling them how to reconstruct one file given the other.
If we can break it down into pieces:

1.  The first line tells us that Mercurial is using the Unix `diff` command
    to compare the last committed and new versions of the file.
2.  The next two lines show us the time stamps of the 2 versions of the file
    that are being compared.
3.  The remaining lines show us the actual differences
    and the lines on which they occur.
    In particular,
    the `+` markers in the first column show where we are adding lines.

Let's commit our change:

~~~
$ hg commit -m "Concerns about Mars's moons on my furry friend"
~~~
{:class="in"}

Checking our project's status:

~~~
$ hg status
~~~
{:class="in"}

we get no output because all of the changes have been committed.
We can see our commits with `hg log`:

~~~
$ hg log
~~~
{:class="in"}
~~~
changeset:   1:9b3b65e50b8c
tag:         tip
user:        Vlad Dracula <vlad@tran.sylvan.ia>
date:        Mon Apr 14 15:52:43 2014 -0400
summary:     Concerns about Mars's moons on my furry friend

changeset:   0:72ab25fa99a1
user:        Vlad Dracula <vlad@tran.sylvan.ia>
date:        Mon Apr 14 14:41:58 2014 -0400
summary:     Starting to think about Mars

~~~
{:class="out"}

Of course sometimes we may not want to commit everything at once.
For example,
suppose we're adding a few citations to our supervisor's work
to our thesis.
We might want to commit those additions,
and the corresponding addition to the bibliography,
but *not* commit the work we're doing on the conclusions
(which we haven't finished yet).
To handle that,
simply do two
(or more)
separate commits,
listing the names of the files to be included in each commit in the `hg commit`
command:

~~~
$ hg commit -m "Cite Frankenstein(2010) and Frankenstein, etal(2011)." methods.txt biblio.txt
...
<later>
...
$ hg commit conclusions.txt -m "Update conclusions re: sunlight."
~~~
{:class="in"}

Notice that the list of file names can come before or after the commit comment
in the `hg commit` command.

Let's add another line to the file for practice and to make our revision
history more interesting:

~~~
$ nano mars.txt
$ cat mars.txt
~~~
{:class="in"}
~~~
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~
{:class="out"}
~~~
$ hg diff
~~~
{:class="in"}
~~~
diff -r 9b3b65e50b8c mars.txt
--- a/mars.txt  Mon Apr 14 15:52:43 2014 -0400
+++ b/mars.txt  Mon Apr 14 16:33:57 2014 -0400
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~
{:class="out"}

So far, so good:
we've added one line to the end of the file
(shown with a `+` in the first column).
Now, let's commit our changes:

~~~
$ hg commit mars.txt -m "Thoughts about the climate"
~~~
{:class="in"}

and look at the history of what we've done so far:

~~~
$ hg log
~~~
{:class="in"}
~~~
changeset:   2:43da31fb96ec
tag:         tip
user:        Vlad Dracula <vlad@tran.sylvan.ia>
date:        Mon Apr 14 16:37:12 2014 -0400
summary:     Thoughts about the climate

changeset:   1:9b3b65e50b8c
user:        Vlad Dracula <vlad@tran.sylvan.ia>
date:        Mon Apr 14 15:52:43 2014 -0400
summary:     Concerns about Mars's moons on my furry friend

changeset:   0:72ab25fa99a1
user:        Vlad Dracula <vlad@tran.sylvan.ia>
date:        Mon Apr 14 14:41:58 2014 -0400
summary:     Starting to think about Mars

~~~
{:class="out"}

#### Exploring History

If we want to see what we changed when,
we use `hg diff` again,
but refer to old versions
using the `--rev` or `-r` flag and the revision numbers:

~~~
$ hg diff --rev 1:2 mars.txt
~~~
{:class="in"}
~~~
diff -r 9b3b65e50b8c mars.txt
--- a/mars.txt  Mon Apr 14 15:52:43 2014 -0400
+++ b/mars.txt  Mon Apr 14 16:44:06 2014 -0400
@@ -1,2 +1,3 @@
 Cold and dry, but everything is my favorite color
 The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~
{:class="out"}
~~~
$ hg diff -r 0:2 mars.txt
~~~
{:class="in"}
~~~
diff -r 72ab25fa99a1 -r 43da31fb96ec mars.txt
--- a/mars.txt  Mon Apr 14 14:41:58 2014 -0400
+++ b/mars.txt  Mon Apr 14 16:37:12 2014 -0400
@@ -1,1 +1,3 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
+But the Mummy will appreciate the lack of humidity
~~~
{:class="out"}

In this way,
we build up a chain of revisions.
The most recent end of the chain is the changeset with the highest revision
number.

To see what changes were made between a particular changeset and its parent
use the `--change` or `-c` flag:

~~~
hg diff --change 1
~~~
{:class="in"}
~~~
diff -r 72ab25fa99a1 -r 9b3b65e50b8c mars.txt
--- a/mars.txt  Mon Apr 14 14:41:58 2014 -0400
+++ b/mars.txt  Mon Apr 14 15:52:43 2014 -0400
@@ -1,1 +1,2 @@
 Cold and dry, but everything is my favorite color
+The two moons may be a problem for Wolfman
~~~
{:class="out"}

#### Recovering Old Versions

All right:
we can save changes to files and see what we've changed---how
can we restore older versions of things?
Let's suppose we accidentally overwrite our file:

~~~
$ nano mars.txt
$ cat mars.txt
~~~
{:class="in"}
~~~
We will need to manufacture our own oxygen
~~~
{:class="out"}

`hg status` now tells us that the file has been changed,
but those changes haven't been committed:

~~~
$ hg status
~~~
{:class="in"}
~~~
M mars.txt
~~~
{:class="out"}

We can put things back the way they were
by using `hg revert`:

~~~
$ hg revert mars.txt
$ cat mars.txt
~~~
{:class="in"}
~~~
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~
{:class="out"}

As you might guess from its name,
`hg revert` reverts to (i.e., restores) an old version of a file.
In this case,
we're telling Mercurial that we want to recover the last committed version
of the file.
If we want to go back even further,
we can use the `--rev` or `-r` flag and a revision number instead:

~~~
$ hg revert --rev 0 mars.txt
~~~
{:class="in"}

Mercurial really doesn't want to cause us to lose our work,
so it defaults to making a backup when we use `hg revert`:

~~~
$ hg status
~~~
{:class="in"}
~~~
? mars.txt.orig
~~~
{:class="out"}

The `mars.txt.orig` file is a copy of `mars.txt` as it stood before the `hg revert` command.
It's not tracked by Mercurial.
It's just there in case we made a mistake and really didn't want to revert,
or in case there's some content from before the revert that we decide that we really do want to copy into `mars.txt`.
When we're sure that we don't need `*.orig` files we can just go ahead and delete them.
If we really don't want Mercurial to create `*.orig` files when we use `hg revert`,
we can use the `--no-backup` option, its short version `-C` or add

~~~
[alias]
revert = revert --no-backup
~~~

to your `~/.hgrc` (or`%USERPROFILE%\Mercurial.ini` if you are using Windows).

The fact that files can be reverted one by one
tends to change the way people organize their work.
If everything is in one large document,
it's hard (but not impossible) to undo changes to the introduction
without also undoing changes made later to the conclusion.
If the introduction and conclusion are stored in separate files,
on the other hand,
moving backward and forward in time becomes much easier.

#### Ignoring Things

What if we have files that we do not want Mercurial to track for us,
like backup files created by our editor
or intermediate files created during data analysis?
Let's create a few dummy files:

~~~
$ mkdir results
$ touch a.dat b.dat c.dat results/a.out results/b.out
~~~
{:class="in"}

and see what Mercurial says:

~~~
$ hg status
~~~
{:class="in"}
~~~
? a.dat
? b.dat
? c.dat
? results/a.out
? results/b.out
~~~
{:class="out"}

Putting these files under version control would be a waste of disk space.
What's worse,
having them all listed could distract us from changes that actually matter,
so let's tell Mercurial to ignore them.

We do this by creating a file in the root directory of our project called `.hgignore`.

~~~
$ nano .hgignore
$ cat .hgignore
~~~
{:class="in"}
~~~
syntax: glob
*.dat
results/
~~~
{:class="out"}

The `syntax: glob` line at the top of the file tells Mercurial that
we want to use the same kind of pattern matching that we use in the shell
(which is known as "globbing" and the patterns as "globs").
The second line tells Mercurial to ignore any file whose name ends in `.dat`
and the third one to ignore everything in the `results` directory.
(If any of these files were already being tracked,
Mercurial would continue to track them.)

Once we have created this file,
the output of `hg status` is much cleaner:

~~~
$ hg status
~~~
{:class="in"}
~~~
? .hgignore
~~~
{:class="out"}

The only thing Mercurial notices now is the newly-created `.hgignore` file.
You might think we wouldn't want to track it,
but everyone we're sharing our repository with will probably want to ignore
the same things that we're ignoring.
Let's add and commit `.hgignore`:

~~~
$ hg add .hgignore
$ hg commit -m "Add the ignore file"
$ hg status
~~~
{:class="in"}

We can also always see the status of ignored files if we want:

~~~
$ hg status --ignored
~~~
{:class="in"}
~~~
I a.dat
I b.dat
I c.dat
I results/a.out
I results/b.out
~~~
{:class="out"}

<div class="keypoints" markdown="1">

#### Key Points
* Create a `~/.hgrc`
  (or `%USERPROFILE%\Mercurial.ini` on Windows) to configure a user name,
  email address, editor, and other preferences once per machine.
* `hg init` initializes a repository.
* `hg status` shows the status of a repository.
* Files are stored in a project's working directory (which users see),
  and the local repository (where committed snapshots are permanently recorded).
* `hg add` tells Mercurial to track files.
* `hg commit` creates a snapshot of the changes to 1 or more files in the local repository.
* Always write a log message when committing changes.
* `hg diff` displays differences between revisions.
* `hg revert` recovers old versions of files.
* The `.hgignore` file tells Mercurial what files to ignore.

</div>

<div class="challenges" markdown="1">

#### Challenges

1.  Create a new Mercurial repository on your computer called `bio`.
    Write a three-line biography for yourself in a file called `me.txt`,
    commit your changes,
    then modify one line and add a fourth and display the differences
    between its updated state and its original state.

</div>
