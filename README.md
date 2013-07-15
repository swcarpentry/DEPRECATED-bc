---
---
Software Carpentry Boot Camps
=============================

This repository's `gh-pages` branch is the starting point for a boot camp website.
You should *not* use this to create the repository students clone for exercises when learning Git:
that should be created from scratch.

Getting Started
---------------

To create a website for a new boot camp:

1.  Create a new repository on GitHub with a name like YYYY-MM-DD-site, e.g., `2014-03-31-ehu`.
    This repository must *not* be a fork of an existing repository.
2.  Clone this repository to your desktop.
3.  Add the template repository `git@github.com:swcarpentry/bc.git` as a remote called `upstream`.
4.  Create a new branch in the desktop clone called `gh-pages`.
5.  Pull content from the GitHub template repository using `git pull upstream gh-pages`.
6.  Edit `index.html` to create the boot camp home page (see below).
7.  Preview your changes (see below).
8.  Push content to your YYYY-MM-DD-site repository using `git push origin gh-pages`.

You may update your boot camp's website whenever you want.

Previewing the Site
-------------------

To preview your boot camp's page(s),
go into its root directory and run:

    jekyll

This will create the directory `./_site` with your rendered pages.

Layout and Variables
--------------------

The `_layouts` directory contains two files:

*   `bootcamp.html`: the layout for bootcamp home pages.
    The material in your bootcamp's `index.html` will be used
    to fill in the {% raw %}{{content}}{% endraw %} section of this page.
*   `lesson.html`: a minimal page for previewing lesson content.

Your bootcamp's `index.html` page (which uses the `bootcamp.html` layout)
must define the following values in its YAML header:

*   `layout` *must* be `bootcamp`.
*   `root` is the path to the repository's root directory (i.e., '.').
    If you create subdirectories for different rooms in your bootcamp,
    or for different lessons,
    set `root` accordingly.
*   `venue` is the name of the institution or group hosting the bootcamp.
*   `address` is the bootcamp venue's street address.
*   `latlng` is the latitude and longitude of the bootcamp site
    (so we can put a pin on our map).
*   `humandate` is the human-friendly dates for the bootcamp (e.g., July 3-4, 2015).
*   `startdate` is the bootcamp's starting date in YYYY-MM-DD format.
*   `enddate` is the bootcamp's ending date in the same format.
*   `registration` is `open` (if anyone is allowed to sign up)
    or `restricted` (if only some people are allowed to take part).
*   `instructor` is a comma-separated list of instructor names.
    This *must* be enclosed in square brackets, as in
    `["Alan Turing", "Grace Hopper"]`.
*   `contact` is the contact email address to use for your bootcamp.

Include Files
-------------

The `_includes` directory contains the following `.html` files:

*   `header.html`: material for the page's head.
*   `banner.html`: the generic banner with the Software Carpentry logo.
*   `footer.html`: the generic footer with links to Software Carpentry's web presence.
*   `javascript.html`: JQuery and Bootstrap Javascript.

The `_includes/bootcamps` directory contains short pieces of standard text
that can be included in boot camp pages using `{% raw %}{% include name.html %}{% endraw %}`:

*   `what.html`: what boot camps are.
*   `who.html`: our intended audience
*   `instructors.html`: creates a list of instructors' names.
*   `python.html`: a brief point-form syllabus for a boot camp using Python.
*   `r.html`: a brief point-form syllabus for a boot camp using R.
*   `requirements.html`: what people need to bring.
*   `contact.html`: how to reach the organizers.

Instructors and Sponsors
------------------------

The `_includes` directory contains directories called `people` and `orgs`
which hold short descriptions of people involved in Software Carpentry
and our financial sponsors respectively.
You can use {% raw %}{% include dir/file.html %}{% endraw %} to include these in your bootcamp home page;
they refer to head shots and logos in `img/people` and `img/orgs`,
and these should be included automatically.
Please send additions to the Software Carpentry administrators.

Lesson Material
---------------

Finally,
the `_includes` directory contains lesson material
that you can incorporate into your pages using `{% raw %}{% include ... %}{% endraw %}`.
Each lesson's material is in its own subdirectory;
different lesson directories are structured differently
but the `db`, `shell`, and `web` directories contain the following files:

*   `instructors.html`: instructors' guide.
*   `opening.html`: opening motivational story.
*   `prereq.html`: discussion of pre-requisites.
*   `reference.html`: a cheat sheet for the subject.
*   `summary.html`: closing summary of the entire lesson.

They also contain subdirectories for various topics,
each of which has the following:

*   `title.md`: the topic title
*   `objectives.html`: the topic's learning objectives
*   `lesson.html`: a long-form prose version of the lesson
*   `summary.html`: the key points of the lesson
*   `challenges.html`: includes all the topic's challenge questions
*   `challenges/some-title.html`: a single challenge question

The images that these files refer to have not yet been created.

Improving This Material
-----------------------

We welcome improvements to the master copy of the boot camp template repository,
particularly new lesson material.
To send them to us:

1.  Fork the `bc` repository on GitHub.
2.  Make that a remote of your desktop copy of your YYYY-MM-DD-site repository.
3.  Push your changes from your desktop to your fork of `bc` on GitHub.
4.  Send a pull request (PR) to the master repository.

FAQ
---

*   *Why does the bootcamp repository have to be created from scratch? Why not fork `bc` on GitHub?*
    <br/>
    Because any particular user can only have one fork of a repository,
    but instructors frequently need to work on several bootcamps at once.

*   *Why are the lesson and topic files HTML instead of Markdown?*
    <br/>
    Primarily convenience---that's what Greg Wilson had in hand to convert.
    These may be converted to Markdown in future.

*   *Why do files like `lessons/db.html` include everything explicitly?*
    <br/>
    Because Liquid does not support parameterized includes like:
    <br/>
    `{% raw %}{% include {{lesson}}/something.html %}{% endraw %}`
    <br/>
    so we can't loop over a set of topics.

*   *Then why use Liquid and Jekyll?  Why not \[some other markup language\] and \[some other converter\]?*
    <br/>
    Because they're the defaults on the site we're encouraging our learners to use.

*   *Where should pages go if multiple boot camps are running at a site simultaneously?*
    <br/>
    Use subdirectories like `2013-07-01-euphoric/beginners`,
    so that main directory names always follow our four-part convention.

*   *How should boot camp instructor lists refer to instructors?*
    <br/>
    Use personal identifiers from the `people` dictionary in `_config.yml`,
    such as `wilson.g`.
    When the boot camp page is expanded by Jekyll on GitHub,
    our template will use this information
    to include names HTML in a standard way.
    A future template will include photos and information about instructors instead.
