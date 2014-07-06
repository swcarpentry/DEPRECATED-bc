Using Git with RStudio
======================

> Git allows groups of people to work on the same documents (often code)
at the same time, and without stepping on each other's toes. It's a
distributed version control system.

(cribbed from [tryGit][trygit])

Intro to practical version control for scientists
=================================================

These [slides][slides] are courtesy of [Bernhard Konrad][konrad].

[konrad]: https://github.com/BernhardKonrad
[slides]: http://htmlpreview.github.io/?https://github.com/BernhardKonrad/2014-02-22-SFU/blob/gh-pages/BK-slides/git-intro.slides.html

Installation and configuration of [git][git]
============================================

If you've already installed and configured [git][git], skip ahead to
[Learn to use git with RStudio](#learngit).

[git]: http://git-scm.com

Install git
-----------

Mac OS 10.9 Mavericks comes with git installed. To check that git is installed,
open a Terminal and run…

```sh
which git
git --version
```

These commands should display something similar to this:

```sh
➜  which git
/usr/bin/git
➜  git --version
git version 1.9.3
```

For all other operating systems, go to the [Git downloads][gitdownloads] web
site, and click on the appropriate icon for your operating system.

If on a Mac the official Git package gives you any trouble, use the following
instructions to install Git using Homebrew.

[gitdownloads]: http://git-scm.com/downloads

Install git using [Homebrew][homebrew]
--------------------------------------

[Homebrew][homebrew] is the missing package manager for Mac OS X. To install
Homebrew and use `brew` to install Git, run…

```sh
ruby -e "$(curl -fsSL https://raw.github.com/Homebrew/homebrew/go/install)"
brew install git
```

Test that git is installed and working by running…

```sh
which git
git --version
```

[homebrew]: http://brew.sh

Configure git
-------------

Git associates your name and e-mail address with each commit, which helps when
multiple people collaborate on a project. To configure your name and e-mail
address in git, open the Terminal and run…

```sh
git config --global user.name 'Your Name'
git config --global user.email 'your@email.com'
```

On a Mac, configure git to remember your password.

```sh
git config --global credential.helper osxkeychain
```

For more help configuring git, see…

+ [GitHub][setupgit]
+ [UBC STAT 540][stat540]

[setupgit]: https://github.com/jennybc/stat540_2014/blob/master/seminars/seminar92_git.md
[stat540]: https://github.com/jennybc/stat540_2014/blob/master/seminars/seminar92_git.md

<a name="configurerstudio">
Configure RStudio to use git
----------------------------
</a>

+ Open RStudio
+ Click *Tools -> Global Options -> Git/SVN*
+ If *Git executable* shows *(none)*, click *Browse* and select the git
  executable installed on your system
  - On a Mac, this will likely be one of
    - `/usr/bin/git`
    - `/usr/local/bin/git`
    - `/usr/local/git/bin/git`
  - On Windows, `git.exe` will likely be somewhere in `Program Files`
+ Click *OK*

<a name="learngit">
Learn to use git with RStudio
=============================
</a>

Create a new project
--------------------

Or if you prefer, see below for instructions to open an existing project.

+ Open RStudio
+ Create a new project
  - Click *File -> New Project -> New Directory -> Empty Project*
  - Check *Create a git repository for this project*

Open an existing project
------------------------

+ Open an existing project
  - Click *File -> Open Project*

If you already have a tab labeled *Git* next to the tabs *Environment* and
*History*, skip these instructions.

+ Enable git for this project
  - Click *Tools -> Version Control -> Project Setup*
  - Click the dropdown box *Version control system* and select *Git*
  - If you don't have a *Git* option go back to [Configure RStudio](#configurerstudio). Do not pass Go. Do not collect $200

Create and commit a file
------------------------

+ Make your first commit
  - Click the *Git* tab
  - Check *Staged* next to `.gitignore` and `hello.Rproj`
  - Click *Commit*
  - Type a message in *Commit message*
  - Click *Commit*
+ Create a new Rmd file
  - Click *File -> New File -> R Markdown*
  - Edit the file and change the title
  - Save the file
+ Commit the new Rmd file
  - Check *Staged* and click *Commit*

Knit the HTML report
--------------------

+ Knit the Rmd file to generate an HTML report
  - Click *Knit HTML*
+ Commit the generated report
  - Check *Staged* for the *md* and *html* files and the *figures* directory
  - Click *Commit*

Change the plot
---------------

+ Replace the *plot* with *ggplot* or *qplot* and save your changes
+ Commit the change
+ Knit the report
+ Commit all the modified files

Make a change and revert it
---------------------------

+ Make an erroneous change to the file and save it
+ Click *Diff* and then *Revert*
+ The erroneous change has been undone and the previous version restored

Delete a file
-------------

+ Create a new file named `doomed.md`
+ Enter some text and save it
+ Delete this doomed file
  - Under the *Files* tab check the box next to `doomed.md`
  - Click *Delete*
+ Under the *Git* tab, a red `D` appears next to the deleted file
+ Stage the change by clicking the checkbox and commit it

Inspect your work
-----------------

+ Make a few more changes and commits
+ Click *History* under the *Git* tab to review your day's work
+ Git has recorded a complete history of your work
+ In the event of impish gnomes introducing errors into your work, you can
  browse through your history, find the gnome to blame, and restore your
  previous good work. Gnomes be damned.

Use the git command line
------------------------

There are many graphical interfaces for git—RStudio is one—but there
is only one git command line interface, which is the common engine
being used behind the scenes. If your graphical interface ever lets
you down, it's useful to peak under the hood.

+ Click *File -> New File -> Text File*
+ Describe your project in this new file
+ Save this file and name it `README.md`
  - Case matters! Name the file `README.md` and not `readme.md` or any
    other variation
  - Don't be imaginative. Get used to being pedantic. Foster your inner OCD
  - md is the extension of a [Markdown](markdown) file
  - Note the yellow question marks indicating the new file that's not
    being tracked by git
+ Open a shell (also known as a Terminal)
  - Under the *Git* tab, click *More -> Shell*
+ Stage `README.md` using the git command line
  - Run `git add README.md`
  - The yellow question mark changes to a green *A*
  - Checking the *Staged* check box in fact runs `git add`
+ Unstage `README.md`
  - Run `git reset README.md`
  - The green *A* changes back to a yellow question mark
  - Unchecking the *Staged* check box in fact runs `git reset`
+ Stage and commit `README.md`
  - Run…

    ```sh
    git add README.md
    git commit -m 'Add README.md'
    ```
  - The `-m` option of `git commit` specifies the git log message
+ Browse the git history in RStudio, and inspect this commit

[markdown]: https://help.github.com/articles/markdown-basics

Learn more about the git command line
-------------------------------------

Go to [tryGit][tryGit] and learn more about the git command line!

[tryGit]: http://try.github.io
