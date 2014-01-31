---
layout: lesson
root: ../../..
title: Git and GitHub
---
**Written by [Matt Davis](mailto:jiffyclub@gmail.com)**

# Git/GitHub

The goal of this lesson is to introduce the students to [Git][git] via
collaboration on [GitHub][].

## Introduction

- Say some introductory stuff about version control in general, and Git/GitHub
  in particular.

*Note: The figures in the [Visual Git Reference][visual git] can be a good
stand-in if you have nowhere to draw.*

## Setup and Signup

- Have everyone configure Git:

        $ git config --global user.name "User Name"
        $ git config --global user.email "user@email.com"
        $ git config --global core.editor "nano"
        $ git config --global color.ui "auto"

- Give a little tour of [GitHub][].
- Have everyone make [GitHub][] accounts.

### Make and Clone

- Make a new demo repo on [GitHub][] explaining the process as you go
  (probably on your personal account).
    - Have [GitHub][] put in a README so it can be cloned.
- Explain that much like a browser navigates to websites using a URL, Git talks
  to remote repositories using a URL.
- Explain the different URL options:
    - Read/write `ssh` protocol requires [ssh keys][], which make it so you
      don't have to enter username/password.
    - `https` protocol takes username/password.
    - `git` protocol is read-only.
- Now we want to get a copy of this down on all of our computers -- `git clone`!
    - Have everyone do this via the `https` URL.
- `ls` to show the new directory and `cd` into it.
- Compare the local listing to what you see on [GitHub][]. Same for now, but
  not for long!

## Basics

### Local Basics

**IMPORTANT:** Make sure you tell people *not* to make their own local changes,
that will make things really complicated later when people pull. Alternatively,
you can go ahead and let them do whatever they like and use it as a teaching
moment on `git reset --hard` in a few minutes when it's time to start the
collaboration.

- On the white board draw a box representing the working area and explain that
  this is where you work and make changes.
- Make a new file called `top-things-to-learn.md` and put the title
  "Top Things We Want to Learn" in it.
- `git status` -- highlight the "Untracked files" section and that Git tells
  you how to add a file to the next commit.

### Composing the Snapshot

- On the white board draw a box representing the staging area (index) and
  explain that this is where we set up the next snapshot of our project.
    - Like a photographer in a studio, we're putting together a shot before
      we actually snap the picture.
    - Connect the working area box and the staging box with `git add`.
- Run `git add top-things-to-learn.md`.
- `git status` -- highlight the "Changes to be committed" section
  and Git telling you how to unstage the new file.

### Taking the Snapshot

- On the white board draw a box representing the project history. Once we take
  a snapshot of the project that snapshot becomes a permanent reference point
  in the project's history that we can always go back to.
    - The history is like a photo album of changes, and each snapshot has a
      time stamp, the name of the photographer, and a description.
    - Connect the staging area to the history with `git commit`.
- Run `git commit` and explain log messages.
    - Summary message at the top, longer one below.
- `git status` -- nothing going on!

### Looking at the History

- `git log` -- Shows description of what we've done.
    - `git log --oneline` -- Abbreviated version.
- Explain the commit hash.
    - Unique identifier of that commit if we ever want to refer to it.
    - Comes from "hashing" stuff associated with the commit, like the changes
      we made.
    - Can demo hashing with Python's `hashlib.sha1`.

### Previewing Changes

- The file we're making is going to be a list of the top things everyone wants
  to learn in the bootcamp. Add your item (e.g. everyone's names) and save.
- `git status` -- point out that now we have a modified file instead of an
  untracked file, but the procedure for adding it to the next snapshot is
  the same.
- Want to see the changes you're about to add? `git diff`!
- `git add`
- `git diff` -- now it doesn't show anything. `git diff` shows differences
  between the working area and the staging area.
    - To see differences between the staging area and the most recent commit
      use `git diff --cached`.
- `git commit -m` -- This time use the `-m` option and show that for short
  commit messages you can just enter them at the command line.

### Undoing Changes

- Say I want to redo the previous commit...
- `git log --oneline` -- grab the commit has for the point we want to go back to.
- `git reset COMMIT`
- `git log --oneline` -- highlight that the latest commit is gone
- `git status` -- But the changes haven't gone anywhere.
- I can now edit the file to fix whatever was wrong and re-commit.

## Sharing

- Now I want to share my changes with everyone so they can start working on
  it too.

### Remotes

- As we said back at the beginning, Git uses URLs to point repositories on other
  computers, in this case [GitHub's][GitHub] servers.
- We can give these remote repositories names so that we don't have to type
  in the full URL all the time, and in fact Git has already set one up for us.
- `git remote` -- there's a remote called "origin".
- `git remote -v` -- we can see that it points to the same URL we cloned from,
  Git sets this up automatically.

### Branches

- On the [GitHub][] view of the repo highlight the branches pull-down -- right
  now it only has one branch called "master", this is another thing Git makes
  for us by default.
- What branch are we on locally? `git branch`.
- Give a short explanation of branches and mention that we will come back to
  them later.
    - Isolated development environments.
- When Git communicates with a remote repository it needs to know what branch
  is changing, in this case we'll just stick with "master".

### Pushing

- Use `push` command to send data to a remote repository, and we also have to
  specify the remote name and branch name: `git push origin master`.
- Refresh the [GitHub][] view.

### Pulling

**IMPORTANT:** If students have been making local commits, this is the time at
which they will need to use `git reset --hard` to get back in sync with you.

- `pull` is the reciprocal command, must specify remote and branch.
- Have everyone `git pull origin master`.

### Collaborate

- Pick a student to add their top thing to learn to the list:
    - Add them to the collaborator list on the demo repo.
    - edit, save, `add`, `commit`, `push`
- Have everyone `pull`.

### Rebase

#### No Conflict

- Have another student add their thing and push.
- Make a change to the README file before pulling.
- Try to push.
- On the white board draw the situation: my repo and the remote repo have
  different development histories and Git doesn't know how to pull things
  together.
- It would be nice if I could move my local change after the remote change.
  (Draw picture.) There's a command for that!
- `git fetch origin` -- This gets information from the remote repository
  but doesn't integrate it with your local repo like `pull` does.
- `git rebase origin/master` -- `origin/master` is how we specify the fetched
  data for the remote named "origin" and it's branch named "master".
    - This replays my local changes on top of the state of the remote repo.
- `git log --oneline` -- Now my change is after the remote change.
- `git push origin master`
- Have everyone pull.

#### With Conflicts

- Again, have a student add their thing and push.
- Before pulling make a change in the same place in the same file.
- Try to rebase as above.
- Explain the conflict message Git prints out.
- Show the conflict messages in the file and how to clean it up.
- Continue the rebase and push the result.
- Have everyone pull.

## Developing in Branches

Often you want to leave a stable version of your code alone while you make some
potentially disruptive changes. Or you and someone else are working on the code
and you want to be able to work without worrying what others are doing.

- It's going to take a long time to get everyone's top thing to learn onto the
  list one at a time, so the room is going to break out into groups and each
  come up with their own list.
- So that they can all do this and then push their result to [GitHub][] each
  is going to work in their own, isolated branch.

### Making a New Branch

*Note: The [Learn Git Branching][] app can be a useful way to
illustrate this part of the lesson.*

- Make a new branch: `git branch matt-list`.
- `git branch` -- highlight the asterisk showing the branch we're currently on.
- Move to that branch: `git checkout matt-list`.
- `git branch` -- asterisk moved!
- Make a change and push.
    - **IMPORTANT:** have to specify new branch named when pushing, not "master".
- `git checkout master` -- show that your change is *not* in master.
- Show how to browse to the other branch on [GitHub][].
- Have each group pick a unique branch name, switch to that branch, and add
  all of their top things to learn to the list and push the result.

### Resolving the Changes

- Browse all of the new branches on [GitHub][].
- Illustrate the situation on the [Learn Git Branching][] app.
- Could just pick one of these branches as the new one "master" and move on,
  but we're more adventurous than that.
- Make sure you're in "master".
- `git fetch origin` -- without a branch name it grabs all of the new branches.
- Pick a branch and `git merge branch-name`.
    - Should be a smooth fast-forward.
    - Illustrate on the [Learn Git Branching][] app.
- Pick another branch and try to merge.
    - Resolve conflicts, add, and commit.
    - Illustrate on the [Learn Git Branching][] app.
- Repeat as necessary.
- Push the final result to [GitHub][].

[git]: http://git-scm.com/
[GitHub]: http://github.com
[ssh keys]: https://help.github.com/articles/generating-ssh-keys
[visual git]: http://marklodato.github.com/visual-git-guide/index-en.html
[Learn Git Branching]: http://pcottle.github.com/learnGitBranching/?NODEMO
