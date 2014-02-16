---
layout: lesson
root: .
title: Setup Instructions
---
These instructions are intended to help learners set up their computers.
Instructors should go through them *before* class starts.

*   [Overview](#overview)
*   [Mac OS X](#macosx)
*   [Windows](#windows)
*   [Linux](#linux)

### Overview
<a name="overview"></a>

#### The Bash Shell

Bash is a commonly-used shell. Using a shell gives you more power to
do more tasks more quickly with your computer.

#### Editor

When you're writing code, it's nice to have a text editor that is
optimized for writing code, with features like automatic color-coding
of key words.

#### Git

Git is a state-of-the-art version control system. It
lets you track who made changes to what when and has
options for easily updating a shared or public version of
your code on [github.com](http://github.com).

#### Python

Python is becoming more and more popular in scientific computing, and
it's a great language for teaching general programming concepts due to
its easy-to-read syntax.  We will be using Python version 2.7.
Installing all the scientific packages for Python individually can be
a bit difficult, so we recommend using an all-in-one installer.

### Mac OS X
<a name="macosx"></a>

#### Bash

The default shell in all versions of Mac OS X is bash, so no need to
install anything.  You access bash from the Terminal (found in
`/Applications/Utilities`).  You may want to keep Terminal in your
dock for this workshop.

#### Editor

We recommend
[Text Wrangler](http://www.barebones.com/products/textwrangler/)
or [Sublime Text](http://www.sublimetext.com/">Sublime Text).

#### Git

Installing Git may require you to first install XCode.
This is a very large download (several gigabytes),
so please do it **before** arriving at the bootcamp.

**For Mac OS X 10.7 and higher:** Go to the
[Xcode website](https://developer.apple.com/xcode/).  Get XCode from
the App Store making certain to install the command line tools (from
the Download preferences pane). Git is included in the command line
tools.

**For Mac OS X 10.6:** If you have Mac OS X 10.6, first get XCode by
going to the
[Apple developer site](https://developer.apple.com/downloads/).
You have to sign in with an Apple ID linked to a Developer account.
If you don't have one, you can register and create one.  Once you log
in, go to page 8 and find "XCode 3.2.6 and iOS SDK 4.3 for Snow
Leopard".  Click to open that section, and then download the `.dmg`
file.  Finally,
[install just git](http://code.google.com/p/git-osx-installer/downloads/list?can=3).

#### Python

We recommend the all-in-one scientific Python installer
[Anaconda](http://continuum.io/downloads.html).
(Installation requires using the shell and if you aren't comfortable
doing the installation yourself just download the installer and we'll
help you at the boot camp.)

*   Download the installer that matches your operating system
    and save it in your home folder.
*   Open a terminal window.
*   Type `bash Anaconda-` and then press tab.
    The name of the file you just downloaded should appear.
*   Press enter.
    You will follow the text-only prompts.
    When there is a colon at the bottom of the screen,
    press the down arrow to move down through the text.
    Type `yes` and press enter to approve the license.
    Press enter to approve the default location for the files.
    Type `yes` and press enter to prepend Anaconda to your `PATH`
    (this makes the Anaconda distribution the default Python).

### Windows
<a name="windows"></a>

#### Software Carpentry Installer

For an all-in-one installer: 

*   Download the [installer](https://raw.github.com/swcarpentry/bc/master/setup/swc-windows-installer.py).
*   If the file opens directly in the browser select **File&rarr;Save Page As** to download it to your computer.
*   Double click on the file to run it.

#### Git Bash

Install Git for Windows by downloading and running
[the installer](https://msysgit.googlecode.com/files/Git-1.8.4-preview20130916.exe).
This will provide you with both Git as well as Bash in the Git Bash program.

#### Editor
  
[Notepad++](http://notepad-plus-plus.org/) is a popular free code editor for Windows.

#### Anaconda Python

*   Download and install [Anaconda CE](http://continuum.io/anacondace.html).
*   Use all of the defaults for installation
    *except* before pressing **Finish**
    make sure to check **Make Anaconda the default Python**.
	
### Linux
<a name="linux"></a>

#### Bash

The default shell is usually `bash`, but if your machine is set up
differently you can run it by opening a terminal and typing `bash`.
There is no need to install anything.

#### Git

If Git is not already available on your machine you can try to install
it via your distro's package manager (e.g. `apt-get`).

#### Editor
  
[Kate](http://kate-editor.org/) is one option for Linux users.

#### Python

We recommend the all-in-one scientific Python installer
[Anaconda](http://continuum.io/downloads.html).
Installation requires using the shell and if you aren't comfortable
doing the installation yourself just download the installer and we'll
help you at the boot camp.
  
*   Download the installer that matches your operating system
    and save it in your home folder.
*   Open a terminal window.
*   Type `bash Anaconda-` and then press tab.
    The name of the file you just downloaded should appear.
*   Press enter. You will follow the text-only prompts. When there is a
    colon at the bottom of the screen press the down arrow to move
    down through the text. Type `yes` and press enter to approve the
    license. Press enter to approve the default location for the
    files. Type `yes` and press enter to prepend Anaconda to your
    `PATH` (this makes the Anaconda distribution the default Python).
