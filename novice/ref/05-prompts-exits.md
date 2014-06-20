---
layout: lesson
root: ../..
title: Recognising prompts and how to exit
---

#### I don't recognise my prompt...where am I?

If the top of your shell window shows...

    GNU nano

...you are in the `nano` text editor.

If your shell window shows...

    ~
    ~
    ~    
    ~
    "filename.txt" ...

...you are in the `vi` text editor.

If the bottom of your shell window shows...

    filename.txt  (Fundamental) ----

...you are in the `emacs` or `xemacs` text editor.

If the bottom of your shell window shows...

    XEmacs: filename.txt  (Fundamental) ----

...you are in the `xemacs` text editor.

If your shell prompt is...

    >>>

... you are in `python`.

If your shell prompt is...

    In [123]:

...you are in `ipython`.

If your shell prompt is...

    >

...you may have typed `'` or `"`, to specify a string, as part of a shell command but have not typed another `'` or `"` to close the string.

If the bottom-left of your shell window shows...

    --More--(...%)

...you are viewing a file using `more`.

If the bottom-left of your shell window shows...

    filename.txt

...you may be viewing a file using `less`.

If the bottom-left of your shell window shows...

    (END)

...you may be viewing a file using `less`.

If the bottom-left of your shell window shows...

    :

...you may be viewing a file using `less` or viewing a `man` page.

If the bottom-left of your shell window shows...

    Manual page...

...you are viewing  a `man` page.

#### How do I exit from...

`nano`:

* Press `CTRL-X`
* If you have unsaved changes, you will be asked to save these - press `y` to save, or `n` to quit without saving.

`vi`:

* Type `:q!` to exit without saving.
* If this text just appears on screen then press `ESC` then type `:q!`

`emacs` or `xemacs`:

* Press `CTRL-X CTRL-C`
* If you have unsaved changes, you will be asked to save these - press `y` to save, or `n` then type `yes` to quit without saving.

`python`:

* Type `exit()` or `CTRL-D`

`ipython`:

* Type `exit()`, or `CTRL-D` then press `y`

`man` page:

* Press `q`

`more` or `less`:

* Press `q`

An open string in a shell, denoted by a `>` prompt:

* If you can see on screen what the character you used to open the strin g was (`'` or `""`) then type the same character again to close the string.
* Or, press `CTRL-C`
