# Undoing changes in Git

## Introduction

Git provides several methods that allow the user to undo changes to files in the working tree and index, and to undo previous commits.
We will cover the techniques based on the `git reset` command.
Note that if all you need to do is go back to a previous commit without undoing changes, `git checkout <commit>` is the command to use.
There will be no actual removal of project history with this method.
For this section we will assume that there is a problem that needs to be undone.

The basic structure of this  command is `git reset <mode> <commit>`.
The `<mode>` option is one of `--soft`, `--mixed`, or `--hard` (omitted defaults to `--mixed`) and tells git to what extent to change the working tree and index.
The `<commit>` option tells git to move the current branch head (`HEAD`) to <commit>.
This has the effect of removing commits from the project history (they will not be seen as output from `git log`).
The table below summarizes the impact of each `<mode>` option.

| Command             | Effect                                         |
|---------------------|------------------------------------------------|
| `git reset --soft`  | No changes to working tree or index            |
| `git reset --mixed` | Remove staged changes from index               |
| `git reset --hard`  | Remove changed files in index and working tree |

Now let's go through an example for each case.

## Undoing changes in the index

Let's say we make some changes to the mars.txt file

```
$ nano mars.txt
$ cat mars.txt

Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
Smooth terrain makes for an easy landing
```

where we have added the last line.
Now we add the new file to the index with the `git add mars.txt` command.

But then we realize that Mars actually has massive volcanoes and does *not* have a smooth surface.
So our mars.txt file has a problem and we do not want to commit the file in this form.
To remove these changes from the index we use the command `git reset --mixed HEAD`.
This specifically tells git to return the index to the state of the last commit, thereby unstaging the bad changes we made to the file.

## Undoing changes in the index and working tree

The previous example removed the bad file from the index but made no changes to the working tree.
With the single command `git reset --hard HEAD` we not only remove the staged changes from the index but also reset the working tree to the state of the last commit.
We have to be careful with this command because one of its important consequences is that the changes in the working tree are lost.

## Undoing bad commits

And finally, we also might want to undo a commit we previously made.
This is especially useful when we realize that the most recent commit needs to include an additional change to the code that we forgot (although in principle we can undo back as many commits as there are in the project history).

Let's look at our project history with the `git log` command.
I will pass the `--oneline` option which just gives the abbreviated SHA-1 number and the commit message, and the `--decorate` option which shows us to which commit the `HEAD` reference points.


```
$ git log --oneline --decorate
005937f (HEAD) Thoughts about the climate
34961b1 Concerns about Mars's moons on my furry friend
f22b25e Starting to think about Mars
```

If we decide that we do not like the changes we have made with the most recent commit and want to reset back to the previous one, we can do so as follows...

```
$ git reset --hard HEAD~
$ git log --oneline --decorate
34961b1 (HEAD) Concerns about Mars's moons on my furry friend
f22b25e Starting to think about Mars
```

The `<commit>` value of `HEAD~` refers to the parent of the last commit.
So this action effectively removed the last commit from the project history and returned the files to the state at the preceding commit.
The `<commit>` option can be set to any value in the project history though all the subsequent work is effectively lost.
Note also that since we selected the `--hard` mode option the index and working tree are also set back to their previous values with that commit.
If we had chosen the `--soft` option, neither the index nor the working tree would have been modified.
