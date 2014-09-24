Software Carpentry Bootcamps
============================

The `bc` repository is the starting point for creating your own bootcamp website:
it contains a template for your bootcamp's home page,
and also hosts the shared lesson materials we have developed.

In this file you will find quick instructions to deploy the homepage of your
bootcamp. If you prefer watch a screencast we have one:
[Setting Up a Software Carpentry Bootcamp Repository](https://vimeo.com/87241285)
(**this can be out of date**). **For more information, please read [README-LONG.md].**

To create a website for a new bootcamp:

1.  Create a [new repository on GitHub](https://github.com/new)
    with a name like YYYY-MM-DD-site, e.g., `2014-03-31-esu`.
    This repository must *not* be a fork of an existing repository.
    Please use the same ID for your bootcamp
    that the Software Carpentry admins are using to identify it
    (i.e.,
    if the admins called the bootcamp `2014-03-31-esu`,
    please *don't* call your repo `euphoric-march-2014`),
    and please use all lower-case
    (i.e., '2014-03-31-esu' instead of '2014-03-31-ESU').

2.  Clone this new repository to your local machine and `cd` into it.
    You can ignore the warning about cloning an empty repository:
    it won't stay empty long.

![Step 1](img/readme/step1.png)

3.  Add the repository `https://github.com/swcarpentry/bc.git` as a remote named `swc`:

    ~~~
    $ git remote add swc https://github.com/swcarpentry/bc.git
    ~~~

![Step 2](img/readme/step2.png)

4.  Create a new branch in the local clone named `gh-pages`.

    ~~~
    $ git checkout -b gh-pages
    ~~~

5.  Pull content from the template repository's `gh-pages` branch into your desktop repository:

    ~~~
    $ git pull swc gh-pages
    ~~~

    This may take a minute or two.

6.  Remove the `swc` remote so that you don't accidentally try
    to push your changes to the main `bc` repository:

    ~~~
    $ git remote rm swc
    ~~~

7.  Edit the lines between `---` at `index.html` with the correct information of
    your bootcamp.

8.  Remove the lines indicated at `index.html`. To know where the block need to
    be delete you can use `grep`:

    ~~~
    $ grep -n 'Remove the block below' index.html
    ~~~

9.  Check if the previous two steps was made correctly. We have a script that
    can do it for you:

    ~~~
    $ make check
    ~~~

9.  Fix the schedule of your bootcamp by editing `_include/schedule.html`.

10. Push content to your YYYY-MM-DD-site repository:

    ~~~
    $ git push origin gh-pages
    ~~~

As soon as your repo has been pushed to GitHub, GitHub will render your pages
at the url:

~~~
http://{your-github-username}.github.io/YYYY-MM-DD-site/
~~~

To update your bootcamp's website just change the files and push again.

Tips
----

1.  The files with the instructions to install all the softwares are at
    `_includes/setup`.
