Contributing Boot Camp Material
===============================

Software Carpentry is an open source/open access project, and we
welcome contributions of all kinds.  By contributing, you are agreeing
that Software Carpentry may redistribute your work under
[these licenses][licenses].  Please see [this page][creators] for
a list of contributors to date.

Improving This Material
-----------------------

We welcome improvements to the master copy of the bootcamp template repository,
particularly new lesson material.
To send them to us as a pull request:

1.  Fork the `bc` repository on GitHub.
2.  Make that a remote named "upstream" of your local YYYY-MM-DD-site repository:

        git remote add upstream https://github.com/<me>/bc.git

(replacing 'me' with your GitHub username)

![Alt text](img/readme/step3.png)

3.  Isolate the changes you want to share in a branch and push them to GitHub
    (you should have added `swcarpentry` as a remote in Step 3 of [Getting Started](#getting-started)):

        git fetch swcarpentry
        git checkout -t swcarpentry/gh-pages -b improvements
        git cherry-pick <commits related to improvements on your gh-pages branch>
        git push upstream improvements

4.  Send a pull request to the master repository on GitHub.

If it is easier for you to send them to us some other way,
please mail us at admin@software-carpentry.org.

Workflow
--------

Software Carpentry uses a development workflow similar to that of
[AstroPy][] and many other open source projects. The AstroPy docs have
excellent sections on:

* [Getting started with git][astropy-git]
* [Developer workflow][astropy-workflow]

File Formats
------------

### Text

Text documents should be in [Markdown][] format and compatible
with [Redcarpet][], the engine GitHub uses to render Markdown.

### Slides

The preferred format for slide presentations is still to be determined.

[AstroPy]: http://astropy.org
[astropy-git]: http://astropy.readthedocs.org/en/latest/development/workflow/index.html#getting-started-with-git
[astropy-workflow]: http://astropy.readthedocs.org/en/latest/development/workflow/development_workflow.html
[creators]: http://software-carpentry.org/badges/creator.html
[licenses]: http://software-carpentry.org/license.html
[Markdown]: http://daringfireball.net/projects/markdown/
[Redcarpet]: https://github.com/vmg/redcarpet
