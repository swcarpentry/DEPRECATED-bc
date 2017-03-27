Getting Started with GitHub
===========================

+ Register an [educational GitHub account][githubedu] if you qualify,
  and a [regular account][github] otherwise
  - A GitHub educational account gives you five private repositories
	for free
  - Public repositories are always free
+ Upload a photo of your self to [Gravatar][gravatar] if you're
  comfortable with that

[github]: https://github.com
[githubedu]: https://education.github.com
[gravatar]: http://gravatar.com

GitHub Documentation
====================

GitHub has great documentation on both git and GitHub.

+ [GitHub Help](https://help.github.com/)
+ [Set Up Git][ghsetup]
+ [Create a Repository][ghcreate]
+ [Fetching a Remote][ghclone]

[ghsetup]: https://help.github.com/articles/set-up-git
[ghcreate]: https://help.github.com/articles/create-a-repo
[ghclone]: https://help.github.com/articles/fetching-a-remote

Set up git
==========

Configure your name and email address if you haven't already.

Open a terminal and run…

```sh
git config --global user.name "Your Name Here"
git config --global user.email "your_email@example.com"
```

On a Mac, also run…

```sh
git config --global credential.helper osxkeychain
```

GitHub has [more detailed instructions][ghsetup] on setting up git.

Create a new RStudio project on GitHub
======================================

+ Open [GitHub][github] in your web browser
+ Create a [new repository][githubnew]
  - Enter a name and description
  - Check *Initialize this repository with a README*
  - Click *Create repository*
+ [Copy the repository URL][ghclone] of your project to the clipboard
  by clicking the *Copy to clipboard* icon next to *HTTPS clone URL*
  on the right-hand side of the project page
+ Open RStudio
+ Create a new project
    - Click *File -> New Project -> Version Control -> Git*
    - Paste the repository URL from the clipboard into *Repository URL*
	- Click *Create Project*
+ Stage the files `.gitignore` and `project.Rproj`
+ Commit these files and push them to GitHub
  - Click *Commit*
  - Click *Pull* and it should respond *Already up to date*
  - Enter a log message and click *Commit*
  - Click *Push*
  - Get in the habit of always clicking *Pull* before commit
  - It will prevent conflicts and save you grief
+ Reload this repository in your web browser and look for these two files
+ Edit `README.md` and add an informative title and description
+ Commit `README.md`
+ Push this commit to GitHub (remember to pull first!)
+ Reload this repository in your web browser and appreciate your beautifully
  rendered `README.md`

[githubnew]: https://github.com/new 

Push an existing RStudio project to GitHub
==========================================

+ Open [GitHub][github] in your web browser
+ Create a [new repository][githubnew]
  - Give the repository the same name as your RStudio project
  - Do ***not*** check *Initialize this repository with a README*
  - Click *Create Project*
+ Copy the two lines of code from the box labeled
  *Push an existing repository from the command line*
+ Open an existing project in RStudio
  - Look for the *Git* tab to ensure it's already using git
+ Open a shell by clicking *Tools -> Shell*
+ Paste the two lines of code that you copied into the shell

  ```sh
  git remote add origin https://github.com/USERNAME/PROJECT.git
  git push -u origin master
  ```
+ Reload this repository in your web browser and browse through your
  handiwork
  - Click on individual files to see their content
  - Click on *commits* to see your history of commits
+ Push your project to GitHub
    - Click *More -> Push Branch*
+ Edit `README.md` and commit the change
  - Remember the mantra: pull, commit, push
+ Reload this repository in your web browser and look for your recent
  change
+ In your GitHub web browser, edit the file `README.md`
  - Click on `README.md`
  - Click *Edit* and make a change
  - Click *Commit changes*
+ In RStudio, pull this change from GitHub
  - Click *More -> Pull Branches*
  - Look for your recent change in RStudio

Learn to use git at the command line
====================================

Learning to use git at the command line is a useful skill to get
yourself out of sticky situations involving conflicts and merges.

[tryGit][trygit] is a fantastic interactive tutorial for learning to
use git at the command line.

[trygit]: http://try.github.io/

Further Reading
===============

## [PHD Comics - "FINAL".doc](http://www.phdcomics.com/comics/archive.php?comicid=1531)

!["FINAL".doc](http://www.phdcomics.com/comics/archive/phd101212s.gif)

## [XKCD - Git Commit](http://xkcd.com/1296/)

![Git Commit](http://imgs.xkcd.com/comics/git_commit.png)

+ [Using Version Control with RStudio][rstudiogit]
+ [An introduction to Git/Github](http://kbroman.github.io/github_tutorial/)
  by Karl Broman, aimed at stats / data science types
+ Ram, K 2013. Git can facilitate greater reproducibility and
  increased transparency in science. Source Code for Biology and
  Medicine 2013 8:7. Go to the
  [associated github repo](https://github.com/karthikram/smb_git)
  to get the PDF (link at bottom of README) and to see a great example
  of how someone managed the process of writing a manuscript with
  git(hub).
+ Blog post [Version control for scientific research](http://blogs.biomedcentral.com/bmcblog/2013/02/28/version-control-for-scientific-research/)
  by Karthik Ram and Titus Brown on the BioMed Central blog February 28
+ Blog post [Getting Started With a GitHub Repository](http://chronicle.com/blogs/profhacker/getting-started-with-a-github-repository)
  from ProfHacker on chronicle.com looks helpful
+ Blog post from The Molecular Ecologist on
  [using GitHub with R and RStudio](http://www.molecularecologist.com/2013/11/using-github-with-r-and-rstudio/)

[rstudiogit]: http://www.rstudio.com/ide/docs/version_control/overview
