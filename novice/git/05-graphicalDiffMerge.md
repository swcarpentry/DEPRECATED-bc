---
layout: lesson
root: ../..
title: Diff and Merge with Graphical Tools
---
<div class="objectives" markdown="1">

#### Objectives
*   Configure Git to use a third party graphical diff/merge tool.
*   Demonstrate how to use third party diff/merge tools in the Git workflow.

</div>

### Intro

Git allows you to integrate 3rd party tools for various tasks. You may have already changed Git's default editor to _nano_.<br>
In this lesson, we will configure Git to use a tool called <a href="https://sourcegear.com/diffmerge/">DiffMerge</a> to help us 
view file changes and resolve merge conflicts.

The Git commands we will be using are `git difftool` and `git mergetool`. 
If you try to use these command before configuring them, Git attempts to use some default diff/merge tools:

~~~
$ git mergetool

This message is displayed because 'merge.tool' is not configured.
See 'git mergetool --tool-help' or 'git help config' for more details.
'git mergetool' will now attempt to use one of the following tools:
tortoisemerge emerge vimdiff
No files need merging

$ git difftool

This message is displayed because 'diff.tool' is not configured.
See 'git difftool --tool-help' or 'git help config' for more details.
'git difftool' will now attempt to use one of the following tools:
kompare emerge vimdiff
~~~
{:class="in"}

So, let's go ahead and configure Git to use [DiffMerge](https://sourcegear.com/diffmerge/)!

### Installation
To install _DiffMerge_ for Windows, OS X, or Linux go to [http://https://sourcegear.com/diffmerge/](http://https://sourcegear.com/diffmerge). There are directions for each system.  

### Integration

To configure Git to use _DiffMerge_, we have to update our `.gitconfig`. 
To do this we simply have to enter a few commands into our command line. 
The exact commands depend on which system you are working on and can be found at 

For Mac: [https://sourcegear.com/diffmerge/webhelp/sec__git__mac.html](https://sourcegear.com/diffmerge/webhelp/sec__git__mac.html)

For Linux: [https://sourcegear.com/diffmerge/webhelp/sec__git__linux.html](https://sourcegear.com/diffmerge/webhelp/sec__git__linux.html)

For Windows: [https://sourcegear.com/diffmerge/webhelp/sec__git__windows.html](https://sourcegear.com/diffmerge/webhelp/sec__git__windows.html)

### Viewing Changes

Remember that you can find information on git commands by adding a `--help` to the end of the command. 
To learn more about `git difftool`, we use the following command:

~~~
$ git difftool --help
~~~
{:class="in"}

From our previous lesson, we have the file `mars.txt`.

~~~
Cold and dry, but everything is my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
~~~
{:class="in"}

Let's alter the file and use `git diff` to see the changes.

~~~
$ nano mars.txt
~~~
{:class="in"}

~~~
Cold and dry, but the rocks are my favorite color
The two moons may be a problem for Wolfman
But the Mummy will appreciate the lack of humidity
Dracula and the Mummy will have to go first, since they don't need oxygen
~~~
{:class="in"}

~~~
$ git diff
diff --git a/mars.txt b/mars.txt
index b36abfd..0e383e3 100644
--- a/mars.txt
+++ b/mars.txt
@@ -1,3 +1,4 @@
-Cold and dry, but everything is my favorite color
+Cold and dry, but the rocks are my favorite color
 The two moons may be a problem for Wolfman
 But the Mummy will appreciate the lack of humidity
+Dracula and the Mummy will have to go first, since they don't need oxygen
~~~
{:class="in"}

Now let's use `git difftool` to see how this looks in _DiffMerge_.

~~~
$ git difftool
Viewing: 'mars.txt'
Launch 'diffmerge' [Y/n]:
~~~
{:class="in"}

By default, Git will prompt you to launch. Type `Y` and hit `Enter`. <br>
If you don't want to be prompted each time you use this command, simply add a `-y` to the end of the command: <br>
`git difftool -y`.<br>

<img src="img/git-difftool.png" alt="Using Diffmerge to view file changes" />
<br>
If you run `git difftool` in a repository which contains changes to multiple files, _DiffMerge_ will be
opened for one file at a time. Once you close the window for the first file, then Git will run
_DiffMerge_ for the next changed file.<br>
This is why it is a good idea to pass the `-y` flag, otherwise you will be prompted at the command line each time
Git tries to open _DiffMerge_ on the next file.<br>


### Resolving Conflicts
Create some conflicts by repeating the steps from the [conflicts lesson](03-conflict.html).

Consider the following conflict:

~~~
$ git diff
diff --cc mars.txt
index 8bbfbcd,2267d7b..0000000
--- a/mars.txt
+++ b/mars.txt
@@@ -1,4 -1,4 +1,8 @@@
  Cold and dry, but everything is my favorite color
++<<<<<<< HEAD
 +The two moons may not be a problem for Wolfman
++=======
+ The two moons will be like Redbull for Wolfman
++>>>>>>> conflict
  But the Mummy will appreciate the lack of humidity
  Dracula and the Mummy will have to go first, since they don't need oxygen
~~~
{:class="in"}

Review the [conflicts lesson](03-conflict.html) for help interpreting the diff above.<br>
Ok, let's use `git mergetool` to see how it displays the conflict above.<br>

<img src="img/git-mergetool.png" alt="Using Diffmerge to view merge conflicts" />

Wait a second, this image shows 3 different panes. The diff above shows two changes, how can there
be 3 panes?
This is why graphical merge tools can be helpful!<br>

A brief summary is contained below.
<br>

<table>
<caption align="bottom"><b><i>Left Pane</i></b></caption>
<tr><td>
<img src="img/git-mergetool-local.png" alt="Local changes" />
</td></tr>
</table>

The left pane shows the local (i.e. your repo) changes which could not be automatically merged.<br>
In this example, this is your local `HEAD` commit.
Notice how the line which contains the changes is marked and color-coded. 

<table>
<caption align="bottom"><b><i>Center Pane</i></b></caption>
<tr><td>
<img src="img/git-mergetool-base.png" alt="Unaltered version" />
</td></tr>
</table>

The center pane shows the _BASE_ file. This is what your file looked like _*before*_ the commit which has the conflict,
which would be the `HEAD~1` commit in this example. 

<table>
<caption align="bottom"><b><i>Right Pane</i></b></caption>
<tr><td>
<img src="img/git-mergetool-remote.png" alt="Remote changes" />
</td></tr>
</table>

The right pane shows the remote _(i.e. the changes you tried to merge into your repo)_ changes which could not be
automatically merged.
<p/>
_DiffMerge_ allows you to pull changes from the left or right to your unchanged version in the middle. 
This is very helpful for complex conflicts.

----

### Summary
We've demonstrated how to integrate a 3rd party diff/merge tool with Git. There are lots of tools which can be 
integrated with Git. The primary reason why we use _DiffMerge_ for this lesson is that it installs and has identical
feature sets for Linux, Mac and Windows.

At this point, it's best to let you explore _DiffMerge_ on your own,
Please refer to the the _DiffMerge_ documentation which may be found through the app itself or online at 
[http://www.sourcegear.com/diffmerge/webhelp/](http://www.sourcegear.com/diffmerge/webhelp/).<br>

Also, `git difftool --help` and `git mergetool --help` contain most (if not all) of the information you will need
to start integrating other tools.<br> 
The same documentation is also available online:<br>
[git difftool](http://git-scm.com/docs/git-difftool)<br>
[git mergetool](http://git-scm.com/docs/git-mergetool)

