Contributing New Material
=========================

Software Carpentry is an open source project,
and we welcome contributions of all kinds:
new and improved lessons,
bug reports,
and small fixes to existing material are all useful.

By contributing,
you are agreeing that Software Carpentry may redistribute your work under
[these licenses](LICENSE.md).

Before beginning anything major,
please read the [README.md](README.md) file in this directory,
which explains how this repository is used to create websites for bootcamps.
We also recommend that you look through [open issues](https://github.com/swcarpentry/bc/issues)
and [pull requests](https://github.com/swcarpentry/bc/pulls) to see what other people could use help with,
and then open an issue of your own to get feedback on your ideas and coordinate with other developers.

**Table of Contents**

*   [Working With GitHub](#working-with-github)  
*   [Previewing](#previewing)
*   [Locations and Formats](#locations-and-formats)
*   [Sample Files](#sample-files)
*   [Labels](#labels)
*   [FAQ](#faq)

Working With GitHub
-------------------

1.  Fork the `swcarpentry/bc` repository on GitHub.

2.  Clone that repository to your own machine.

3.  Create a branch from `master` for your changes.
    Give your branch a meaningful name,
    such as `fixing-typos-in-novice-shell-lesson`
    or `adding-tutorial-on-visualization`.

4.  Make your changes, commit them, and push them to your repository on GitHub.

5.  Send a pull request to the `master` branch of `[swcarpentry/bc](http://github.com/swcarpentry/bc)`.

If it is easier for you to send them to us some other way,
please mail us at
[admin@software-carpentry.org](mailto:admin@software-carpentry.org).
Given a choice between you creating content or wrestling with Git,
we'd rather have you doing the former.

Previewing
----------

To preview changes before committing,
run the command `make site`.
This runs Jekyll with the same flags that GitHub uses when things are committed to the `gh-pages` branch
and puts the results in a directory called `_site`.


You should also run `make check` before pushing changes to your `index.html` home page
to your repository.
If you don't have Make installed,
you can run the same checks using:

~~~
python bin/swc_index_validator.py ./index.html
~~~

This checks that the bootcamp's instructors are listed,
that a contact email address has been set up,
and so on.

Other useful commands in the main Makefile are:

*   `make commands` (or just `make` on its own): list available commands.
*   `make clean`: remove editor backup files and the `_site` directory.
*   `make fixme`: list uses of the word `FIXME` in source files.

The commands to convert IPython Notebooks to Markdown
are stored in a separate Makefile called `ipynb.mk`
to simplify maintenance
and ensure that the main Makefile only does what Jekyll on GitHub will do.
To re-do conversion of notebooks to Markdown files,
use `make ipynb`.

Locations and Formats
---------------------

Every subject has a sub-directory of its own,
while individual topics are files in that directory.
For example,
the `novice/git` directory holding our introduction to Git for newcomers
contains the files
`00-intro.md`,
`01-backup.md`,
and so on.
(We use two digits followed by a one-word topic key
to ensure files appear in the right order when listed.)

Lessons may be written in Markdown,
as IPython Notebooks,
or in other formats.
However,
as explained in [the README file](README.md),
Jekyll (the tool GitHub uses to create websites)
only knows how to handle Markdown and HTML.
if some other format is used,
the author of the lesson must
add the generated Markdown to the repository.
This ensures that people who *aren't* familiar with some format
don't have to install the tools needed to work with it
(e.g.,
R programmers don't have to install the IPython Notebook).

> If a lesson is in a format we don't already handle,
> the author must also add something to the Makefile
> to re-create the Markdown from the source.
> Please check with us if you plan to do this.

Sample Files
------------

The directory `misc` contains files that can be used as starting points for lessons.
These files explain how to format titles,
objectives,
exercises,
key points,
and code fragments.
In addition,
the IPython Notebook file has metadata in various cells
to ensure that the generated HTML pages have the right style.
If you are creating a new lesson,
please copy one of these files to use as a starting point.

Labels
------

We use labels to categorize new Issues and Pull Requests. If you are searching
for Issues to help resolve or Pull Requests to review, you can filter with these
labels to help find subjects you are interested in. For example, if you are
interested in improving the novice Python lessons, you can filter with the
labels "Python" and "novice".

If you have "Contributor" status for `swcarpentry/bc`, please help organize
everything by adding informative labels to Issues and Pull Requests you create
as well as Issues and Pull Requests made by others.

See the list below for descriptions of some common labels:

*   bug: something's wrong and needs to be fixed
*   discussion: this issue is being used as a mini-mailing list to discuss something of specialized interest
*   enhancement: identifies improvements to existing material (rather than bug fixes)
*   getting-started: the task is suited for a newcomer, e.g. someone in instructor training
*   intermediate: concerns the intermediate set of teaching material
*   lessons: anything related to lessons that doesn't fit into a more specific category
*   novice: concerns the novice set of teaching material
*   question: asking for help or getting feedback on an idea
*   tools: concerns the tools used to build and manage lesson material
*   work in progress: the pull request is not ready to be merged

FAQ
---

*   *Where can I get help?*
    <br/>
    Mail us at [admin@software-carpentry.org](mailto:admin@software-carpentry.org),
    come chat with us on [our IRC channel](irc://moznet/sciencelab),
    or join our [discussion list](http://software-carpentry.org/pages/discuss.html)
    and ask for help there.
