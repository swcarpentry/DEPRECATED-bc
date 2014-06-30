---
layout: lesson
root: ../..
title: Recognising prompts and how to exit
---

#### I don't recognise my prompt...where am I?

-----

**If the top of your shell window shows...**

    GNU nano

...you are in the `nano` text editor.

**To exit `nano`:**

* Press `CTRL-X`
* If you have unsaved changes, you will be asked to save these - press `y` to save, or `n` to quit without saving.

-----

**If your shell window shows...**

    ~
    ~
    ~    
    ~
    "filename.txt" ...

...you are in the `vi` text editor.

**To exit `vi`:**

* Press `ESC` then type `:q!` to exit without saving.
* To save your changes, after `ESC`, type `:wq`
-----

**If the bottom of your shell window shows...**

    filename.txt  (Fundamental) ----
or:

    XEmacs: filename.txt  (Fundamental) ----

...you are in the `emacs` or `xemacs` text editor.

**To exit `emacs` or `xemacs`:**

* Press `CTRL-X CTRL-C`
* If you have unsaved changes, you will be asked to save these - press `y` to save, or `n` then type `yes` to quit without saving.

-----

**If your shell prompt is...**

    >>>

... you are in `python`.

**To exit `python`:**

* Type `exit()` or `CTRL-D`

-----

**If your shell prompt is...**

    In [123]:

...you are in `ipython`.

**To exit `ipython`**:

* Type `exit()`, or `CTRL-D` then press `y`

------

**If your shell prompt is...**

    >

...you may have typed `'` or `"`, to specify a string, as part of a shell command but have not typed another `'` or `"` to close the string. You may also have an open parenthesis `(`

**To recover from an unterminated expression:** 

* If you can see on screen what the character you used to open the string was (`'` or `""`) then type the same character again to close the string.
* Or, press `CTRL-C`

------
**If the bottom-left of your shell window shows...**

    --More--(...%)
or

    filename.txt
or

    (END)
or

    :

...you may be viewing a file using `more` or `less` or viewing a `man` page.

**To exit** `more` or `less` or a `man` page:

* Press `q`

