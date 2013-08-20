---
layout: lesson
root: ../..
title: Installing a Python Package
---
One of my favorite tools is a Python tool called **grin**.  If you type grin at the command prompt, what do you get?

    No command `grin` found did you mean:
    ...

Looks like it is not installed, but we can fix this.

**grin** is a tool written in Python.
We will want to install it, so we will use a simple package manager for Python to install it called **pip**

    pip install grin

Hmmmm, it doesnt like us, and we get a permission denied error.
We could ask for the sys-admin to install for us,
but we are impatient,
and we now have tools to get this to work.
Here is what we can do:

1. make a directory called `local` in our home direcotry
2. use pip to install `grin` to this local directory using **--install-options="--prefix=/home/swc/local"**
3. add `/home/swc/local/bin` to our `PATH` so we can use grin (the binary **grin** was installed in `/home/swc/local/bin`)
4. add `/home/swc/local/lib/python2.7/site-packages` to our `PYTHONPATH`

Here are the commands

    cd /home/swc
   
    mkdir local
   
    pip install --install-option='--prefix=/home/swc/local' grin

but because this is not defined in **$PATH** the system can *NOT* find it.
So for your last part of the day,
lets look at a text editor.

### gedit

**gedit** is a simple text editor,
much like Microsoft Word is a text editor,
but it is specially designed to work with code and scripts.
To open gedit,
type `gedit` at the command line:

    gedit

This will open a blank text editor.
If you use the menu to open `generate_data.py`,
it will open the text file.
Two important things to notice:

1. The text is highlighted,
   and you will see a `python` tab at the bottom of the window.
   This lets you know that it is higlighting the text in the file as if it were Python code.
   You will find this very helpful in the next session.
2. There will be line numbers on the left side of the page.
   This will also help in communicating to others where you are in the document.
   If I tell you to look at line 20,
   you should see that this line contains the text

    birthmonths= range(1,13)

To close gedit, use the menu `File -> Quit`,
or click on the `X` in the upper right hand corner of the window.

### Hidden files revisited

Remember when we used `ls -a` to find hidden files.
We are going to edit one of those hidden files.
First go to your home directory `/home/swc`
(a quick way to do this is just type `cd` as it always takes you home.
Use `pwd` to verify where you are:

    cd
    pwd

We are going to update your Bash resource file.

### .bashrc

The **.bashrc** file in your home directory controls the behavior of your shell.
Our goal is to update it so the system will look in `/home/swc/local/bin` for executable programs,
and in `/home/swc/local/lib/python2.7/site-packages` for your newly installed Python module.
This is done by adding `/home/swc/local/bin` to your *$PATH* environment variable
and `/home/swc/local/lib/python2.7/site-packages` to your *$PYTHONPATH* variable.

So lets use gedit to open our `.bashrc` file.

    gedit .bashrc

There is no syntax higlighting.
It would be helpful to have some,
so use your mouse to click on **Plain Text**,
then choose the **sh** interpreter to highlight your text.
There is alot of text in this file,
and we do not have time to cover it in this course.
For now we will just scroll to the bottom of the text file and add a couple new lines.

To update your **PATH** and **PYTHONPATH** environment variables, add the following lines:

    export PATH=$PATH:/home/swc/local/bin
    export PYTHONPATH=$PYTHONPATH:/home/swc/local/lib/python2.7/site-packages

NOTE:: 
1. There are no spaces between `PATH` and `=`.
2. There is a colon **:** between `$PATH` and `/home/swc/local/bin`.
3. **$PATH** gets the current value of PATH,
   and makes sure you dont lose these paths.
   If you did not do this, many other programs would stop working,
   as your system would not know where to find them.

Save your edited file `File -> Save`.
For your changes to take effect,
you need to open a new shell.

### grin

Now when you type **grin** at the prompt you should see this:

    usage: grin [-h] [-v] [-i] [-A AFTER_CONTEXT] [-B BEFORE_CONTEXT] [-C CONTEXT]
                [-I INCLUDE] [-n] [-N] [-H] [--without-filename] [--emacs] [-l]
                [-L] [--no-color] [--use-color] [--force-color] [-s]
                [--skip-hidden-files] [-b] [--skip-backup-files] [-S]
                [--skip-hidden-dirs] [-d SKIP_DIRS] [-D] [-e SKIP_EXTS] [-E]
                [--no-follow] [--follow] [-f FILE] [-0] [--sys-path]
                regex [files [files ...]]
    grin: error: too few arguments

Congradulate yourself, you just installed a Python module and it works!

Now lets see why **grin** is so cool.
Basically it lets you search recursively into directories for a pattern in text files.
Go back to your project directory:

    cd boot-camps/shell/ImplantProject

Lets look for the text `prince` in any file in the `rawdata` directory

    grin prince rawdata/

You should see

1. a list of files that have prince in them along with the line number where the text was found, and
2. each instance of `prince` is highlighted.
