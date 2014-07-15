---
layout: lesson
root: ../..
title: Forking a Repository
---

The model shown in the [main lesson](../git/02-collab.html)
in which everyone pushes and pulls from a single repository,
is perfectly usable,
but it will only work if you have write access to the repository.
Sometimes you will want to contribute to someone else's repository
and you won't be able to push your changes to it.
Instead, you can create your own copy of the repository on Github,
push your changes to your copy,
and ask the original author to review and possibly accept your changes
back into the original repository.

Suppose Wolfman wants to be able to make changes to Dracula's project on Github.
Instead of creating a new project,
Wolfman [forks](../../gloss.html#fork) it,
i.e., clones it on GitHub. He does this using the GitHub web interface:

<img src="img/git-fork-ui.png" alt="The Fork Button" />

He then [clones](../../gloss.html#repository-clone) his own GitHub repository,
not Dracula's,
to give himself a desktop copy:

<img src="img/git-forking-01.svg" alt="After Forking on GitHub" />

Now Wolfman can make changes locally, on his machine.
He can make a change to the project
and commit it to his local repository.
Then he can use `git push` to copy those changes to GitHub:

<img src="img/git-forking-02.svg" alt="After Pushing to Fork" />

Once on GitHub
the changes are shared with new potential collaborators.
Wolfman can share his fork
or others can fork his repository with its changes.
He can even share his changes with Dracula.
Dracula can view them
and offer feedback
or accept the changes
and merge them into his own repository,
just as he would review Wolfman's paper before publishing it.

Likewise, Wolfman can review future changes Dracula makes
and decide to merge them into his fork of the project.
If Dracula decides to remove an important feature,
or rewrite the project,
or delete the project entirely,
Wolfman retains his unmodified copy.
Most importantly,
Dracula does not have to give Wolfman access
to make changes to his repository
in order for Wolfman to create
and share
useful modifications.

This review and merge process
can be done manually
by working in branches,
using `git pull` to get changes from one repository,
and `git push` to apply the changes to another repository.
GitHub has a tool for proposing two repositories share changes:
the [pull request](../../gloss.html#pull-request).

In order to share his changes with Dracula,
Wolfman can create a [pull request](../../gloss.html#pull-request),
which notifies Dracula that Wolfman wants to merge some changes into Wolfman's repository:

<img src="img/git-forking-03.svg" alt="After Creating Pull Request" />

A pull request is a merge waiting to happen.
When Dracula views it online,
he can see and comment on the changes Wolfman wants to make.
Wolfman can update his branch on his fork
and the pull request will automatically update with the changes.
Likewise, Dracula can update his branch
and the pull request will update to reflect the changes.
When Dracula likes the changes
and wants to merge them into his project
he can do so with the click of a button:

<img src="img/github-merge-ui.png" alt="Mergeing a Pull Request" />

Dracula may take several days to get around to reviewing Wolfman's Pull Request.
Rather than being stalled until then,
Wolfman can switch to another branch and work on something else,
then switch back to make modifications if necessary.
Once the changes are accepted,
Wolfman can delete the branch; provided it has been merged into `master`
(or some other branch),
the only thing that will be lost is the pointer with the branch name,
not the changes themselves.
