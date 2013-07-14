---
---
Software Carpentry Boot Camps
=============================

This repository's `gh-pages` branch is the default starting point for a bootcamp's web site.
It is also a submodule of the `site` repository that holds Software Carpentry's main web site
so that the site-building script can harvest metadata from pages
for inclusion in the boot camp overview page.

Getting Started
---------------

The admins will normally do the following for you:

*   Clone this repository to create one that you can use for your bootcamp.

Once this is done, you should:

*   Fill in the header of the `index.html` page with your bootcamp's location, dates, instructors, etc.
*   Include relevant setup instructions from `_includes/setup/*.html`.
*   Preview your changes locally.
*   Push your changes to GitHub to update your bootcamp's public web site.
*   Send pull requests back to the master copy whenever you have changes to offer.

Previewing the Site
-------------------

To preview your boot camp's page(s),
go into its root directory and run:

    jekyll

This will create the directory `./_site` with your rendered pages.

Stock Material
--------------

The `_includes` directory contains short pieces of standard text
that can be included in boot camp pages using `{% raw %}{% include name.html %}{% endraw %}`:

*   `what.html`: what boot camps are.
*   `who.html`: our intended audience
*   `instructors.html`: creates a list of instructors' names.
*   `python.html`: a brief point-form syllabus for a boot camp using Python.
*   `r.html`: a brief point-form syllabus for a boot camp using R.
*   `requirements.html`: what people need to bring.
*   `contact.html`: how to reach the organizers.

Lesson Material
---------------

The `_includes/lessons` directory includes standard lesson material.
To incorporate material from it,
you may:

*   copy and paste it into your boot camp's `index.html` page;
    or
*   copy and paste the page into your boot camp's directory,
    then use `{% raw %}{% include ... %}{% endraw %}`
    to include it in your `index.html`.

You should only use the first option if you are modifying the material.

Each lesson's directory contains the following files,
which you can include using `{% raw %}{% include lesson/whatever.html %}{% endraw %}`:

*   `instructors.html`: instructors' guide.
*   `opening.html`: opening motivational story.
*   `prereq.html`: discussion of pre-requisites.
*   `reference.html`: a cheat sheet for the subject.
*   `summary.html`: closing summary of the entire lesson.

You may also use `{% raw %}{% include lesson/topic/whatever.html %}{% endraw %}`
to include topic material in your page(s).
`lesson` is a lesson name, such as `db`;
`topic` is the particular topic, such as `select`;
and `whatever` is one of:

*   `title.md`: the topic title
*   `objectives.html`: the topic's learning objectives
*   `lesson.html`: a long-form prose version of the lesson
*   `summary.html`: the key points of the lesson
*   `challenges.html`: includes all the topic's challenge questions
*   `challenges/some-title.html`: a single challenge question

FAQ
---

*   *Why are the lesson and topic files HTML instead of Markdown?*
    <br/>
    Primarily convenience---that's what Greg Wilson had in hand to convert.
    These may be converted to Markdown in future.

*   *Why does `lessons/db.html` include everything explicitly?*
    <br/>
    Because Liquid does not support parameterized includes like:
    <br/>
    `{% raw %}{% include {{lesson}}/something.html %}{% endraw %}`
    <br/>
    so we can't loop over a set of topics.

*   *Then why use Liquid and Jekyll?  Why not \[some other markup language\] and \[some other converter\]?*
    <br/>
    For the same reason we use Git:
    they're the defaults on the site we're encouraging our learners to use.

*   *Where should pages go if multiple boot camps are running at a site simultaneously?*
    <br/>
    Use sub-directories like `2013-07-01-euphoric/beginners`,
    so that main directory names always follow our four-part convention.

*   *How should boot camp instructor lists refer to instructors?*
    <br/>
    Use personal identifiers from the `people` dictionary in `_config.yml`,
    such as `wilson.g`.
    When the boot camp page is expanded by Jekyll on GitHub,
    our template will use this information
    to include names HTML in a standard way.
    A future template will include photos and information about instructors instead.

*   *Where should new material be put?*
    <br/>
    If you would like to add material to `lessons`,
    or modify material that is there,
    please send a pull request rather than editing directly:
    other people may be depending on its current state.
