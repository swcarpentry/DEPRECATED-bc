---
layout: lesson
root: ../..
title: Collaborating
---
<div class="objectives" markdown="1">

#### Objectives
* Explain what remote repositories are and why they are useful.
* Explain what happens when a remote repository is cloned.
* Explain what happens when changes are pushed to or pulled from a remote repository.

</div>

Version control really comes into its own
when we begin to collaborate with other people.
We already have most of the machinery we need to do this;
the only thing missing is to copy changes from one person's repository to another.

Systems like Mercurial allow us to move work between any two repositories.
In practice,
though,
it's easiest to use one copy as a central hub,
and to keep it on the web rather than on someone's laptop.
Most programmers use hosting services like [BitBucket](http://bitbucket.org)
and similar to hold those master copies;
we'll explore the pros and cons of this in the final section of this lesson.

Let's start by sharing the changes we've made to our current project with the world.
Log in to BitBucket,
then click on the icon in the top right corner to create a new repository called `planets`:

<img src="img/bitbucket-create-repo-01.png" width="600px" alt="Creating a Repository on BitBucket (Step 1)" />

Name your repository "planets" and then click "Create Repository":

<img src="img/bitbucket-create-repo-02.png" width="600px" alt="Creating a Repository on BitBucket (Step 2)" />

As soon as the repository is created,
BitBucket displays a page with a URL and some information on how to configure your local repository:

<img src="img/bitbucket-create-repo-03.png" width="800px" alt="Creating a Repository on BitBucket (Step 3)" />

This effectively does the following on BitBucket's servers:

~~~
$ mkdir planets
$ cd planets
$ hg init
~~~
{:class="in"}

Our local repository still contains our earlier work on `mars.txt`,
but the remote repository on BitBucket doesn't contain any files yet:

<img src="img/hg-freshly-made-bitbucket-repo.svg" alt="Freshly-Made BitBucket Repository" />

The next step is to connect the two repositories.
We do this by making the BitBucket repository a [remote](../../gloss.html#repository-remote)
for the local repository.
The home page of the repository on BitBucket includes
the string we need to identify it after clicking on "I have an existing project to push up":

<img src="img/bitbucket-find-repo-string.png" width="800px" alt="Where to Find Repository URL on BitBucket" />

Change the 'ssh://' string to 'https://' in the url [protocol](../../gloss.html#protocol).
It's slightly less convenient for day-to-day use,
but much less work for beginners to set up.

Copy that URL from the browser,
go into the local `planets` repository,
and use your text editor to add the following lines to `.hg/hgrc`:

~~~
[paths]
default = https://bitbucket.org/vlad/planets
~~~

Make sure to use the URL for your repository rather than Vlad's;
the only difference should be your username instead of `vlad`.

We can check that the command has worked by running `hg paths`:

~~~
$ hg paths
~~~
{:class="in"}
~~~
default = https://bitbucket.org/vlad/planets
~~~
{:class="out"}

Once the default path is set up,
this command will push the changes from our local repository
to the repository on BitBucket:

~~~
$ hg push
~~~
{:class="in"}
~~~
pushing to https://bitbucket.org/vlad/planets
searching for changes
adding changesets
adding manifests
adding file changes
added 1 changesets with 1 changes to 1 files
~~~
{:class="out"}

Our local and remote repositories are now in this state:

<img src="img/bitbucket-repo-after-first-push.svg" alt="BitBucket Repository After First Push" />

We can pull changes from the remote repository to the local one as well:

~~~
$ hg pull
~~~
{:class="in"}
~~~
pulling from https://bitbucket.org/vlad/planets
searching for changes
no changes found
~~~
{:class="out"}

Pulling has no effect in this case
because the two repositories are already synchronized.
If someone else had pushed some changes to the repository on BitBucket,
though,
this command would download them to our local repository.

We can simulate working with a collaborator using another copy of the repository on our local machine.
To do this,
`cd` to the directory `/tmp`.
(Note the absolute path:
don't make `tmp` a subdirectory of the existing repository).
Instead of creating a new repository here with `hg init`,
we will [clone](../../gloss.html#repository-clone) the existing repository from BitBucket:

~~~
$ cd /tmp
$ hg clone https://bitbucket.org/vlad/planets
~~~
{:class="in"}

`hg clone` creates a fresh local copy of a remote repository.
(We did it in `/tmp` or some other directory so that we don't overwrite our existing `planets` directory.)
Our computer now has two copies of the repository:

<img src="img/hg-after-duplicate-clone.svg" alt="After Creating Duplicate Clone of Repository" />

Let's make a change in the copy in `/tmp/planets`:

~~~
$ cd /tmp/planets
$ nano pluto.txt
$ cat pluto.txt
~~~
{:class="in"}
~~~
It is so a planet!
~~~
{:class="out"}
~~~
$ hg add pluto.txt
$ hg commit -m "Some notes about Pluto"
~~~
{:class="in"}

then push the change to BitBucket:

~~~
$ hg push
~~~
{:class="in"}
~~~
pushing to https://bitbucket.org/vlad/planets
searching for changes
adding changesets
adding manifests
adding file changes
added 1 changesets with 1 changes to 1 files
~~~
{:class="out"}

Note that we didn't have to create a remote called `default`:
Mercurial does this automatically,
using that name,
when we clone a repository.
(This is why `default` was a sensible choice earlier
when we were setting up remotes by hand.)

Our three repositories now look like this:

<img src="img/hg-after-change-to-duplicate-repo.svg" alt="After Pushing Change from Duplicate Repository" />

We can now download changes into the original repository on our machine:

~~~
$ cd ~/planets
$ hg pull
~~~
{:class="in"}
~~~
pulling from https://bitbucket.org/vlad/planets
searching for changes
adding changesets
adding manifests
adding file changes
added 1 changesets with 1 changes to 1 files
~~~
{:class="out"}

which gives us this:

<img src="img/hg-after-pulling-to-local-repo.svg" alt="After Pulling Change to Local Repository" />

In practice,
we would probably never have two copies of the same remote repository
on our laptop at once.
Instead,
one of those copies would be on our laptop,
and the other on a lab machine,
or on someone else's computer.
Pushing and pulling changes gives us a reliable way
to share work between different people and machines.

<div class="keypoints" markdown="1">

#### Key Points
* A local Mercurial repository can be connected to one or more remote repositories.
* Use the HTTPS protocol to connect to remote repositories until you have learned how to set up SSH.
* `hg push` copies changes from a local repository to a remote repository.
* `hg pull` copies changes from a remote repository to a local repository.
* `hg clone` copies a remote repository to create a local repository with a remote called `default` automatically set up.

</div>

<div class="challenges" markdown="1">

#### Challenges

1.  Create a repository on BitBucket,
    clone it,
    add a file,
    push those changes to BitBucket,
    and then look at the [timestamp](../../gloss.html#timestamp) of the change on BitBucket.
    How does BitBucket record times, and why?

</div>
