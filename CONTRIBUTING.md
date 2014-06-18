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
read the [README.md](README.md) file in this directory.
It explains how this repository is used to create websites for bootcamps.
We also recommend that you open an issue in the `swcarpentry/bc` Issue Tracker
to get feedback on your ideas and coordinate with other developers.

**Table of Contents**

*   [Working With GitHub](#working-with-github)
*   [Locations](#locations)
*   [Formats](#formats)
*   [Sample Files](#sample-files)
*   [Previewing](#previewing)
*   [FAQ](#faq)

Working With GitHub
-------------------

In general, we are using Git and GitHub for almost everything.
But if it is easier for you to send your contributions to us some other way,
please mail us at
[admin@software-carpentry.org](mailto:admin@software-carpentry.org).
Given a choice between you creating content or wrestling with Git,
we'd rather have you doing the former.

The GitHub workflow looks like this:

1.  Fork the `swcarpentry/bc` repository on GitHub.

2.  Clone that repository to your own machine.

3.  Create a branch from `master` for your changes.
    Give your branch a meaningful name,
    such as `fixing-typos-in-novice-shell-lesson`
    or `adding-tutorial-on-visualization`.

4.  Make your changes, commit them, and push them to your repository on GitHub.

5.  Send a pull request to the `master` branch of `[swcarpentry/bc](http://github.com/swcarpentry/bc)`.

Locations
---------

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

Formats
-------

We are using [Jekyll](https://en.wikipedia.org/wiki/Jekyll_%28software%29)
to convert lessons in [Markdown](https://en.wikipedia.org/wiki/Markdown)
to [HTML](https://en.wikipedia.org/wiki/HTML).
Therefore, we require every lesson to be available in Markdown.
But this does not mean that you have to write it in Markdown as well.

The following list contains more detailed instructions for different formats:

*   Markdown: Nothing to do.
*   IPython Notebook: Generate and provide a Markdown version in addition
    to the IPython Notebook.
    Look at our existing files to get an idea how the IPython Notebook
    should look like.
*   Other formats: Generate and provide a Markdown version in addition
    to the file in the original format.
    In order to be able to recreate the Markdown from source,
    include the conversion from this other format to Markdown
    in the `Makefile` as well. This ensures that people who *aren't* familiar
    with some format don't have to install the tools needed to work with it
    (e.g., R programmers don't have to install the IPython Notebook).
    **Please check with us if you plan to do this.**

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

Previewing
----------

To preview changes before committing,
run the command `make site`.
This runs Jekyll with the same flags that GitHub uses when things are committed to the `gh-pages` branch.
Jekyll's output is stored in a directory called `_site`.

Other useful commands in the main `Makefile` are:

*   `make commands` (or just `make` on its own): list available commands.
*   `make check`: check that the repository's `index.html` file is properly formatted.
*   `make clean`: remove editor backup files and the `_site` directory.
*   `make fixme`: list uses of the word `FIXME` in source files.

The commands to convert IPython Notebooks to Markdown
are stored in a separate `Makefile` called `ipynb.mk`
to simplify maintenance
and ensure that the main `Makefile` only does what Jekyll on GitHub will do.
To re-do conversion of notebooks to Markdown files,
use `make ipynb`.

FAQ
---

*   *Where can I get help?*
    Mail us at [admin@software-carpentry.org](mailto:admin@software-carpentry.org),
    come chat with us on [our IRC channel](irc://moznet/sciencelab),
    or join our [discussion list](http://software-carpentry.org/pages/discuss.html)
    and ask for help there.
*   *Is there a list of open issues to tackle?*
    The [open issues for the bc repository](https://github.com/swcarpentry/bc/issues)
    and also the [related projects](http://www.software-carpentry.org/pages/create.html#related)
    listed on the create page are a good starting point.

