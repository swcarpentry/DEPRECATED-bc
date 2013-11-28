---
layout: lesson
root: ../..
title: Introduction to R
---
This repository gives the R content from Duke software Carpentry May 20 - 21 2013, which was very similar to the NESCent bootcamp which ran May 17 - 18,2013.

Here are pages with set-up instructions, etc.:

*   [NESCent bootcamp](http://swcarpentry.github.io/bootcamps/2013-05-16-nescent/)
*   [Duke bootcamp](http://swcarpentry.github.io/bootcamps/2013-05-20-duke/)

Day 1 was essentially all R, using RStudio. Main instructor Jenny Bryan. At the end of the day, we had some input data (an excerpt from the Gapminder data), several R scripts, and various outputs that these scripts left behind, e.g. figures, numerical results, analytical reports. The first commit to this repo is a snapshot of where we are at the end of Day 1. The README.md associated with the `code` subdirectory explains what each file does.

Our goal was to introduce students to data analysis in R, so focus was not on progamming *per se*. Coverage of writing functions, control structures, package development etc. is essentially non-existent. That was intentional. Had to do some visualization, which we did with `lattice`, since that's what Jenny knows. More graphics would be good -- but how to fit in? Also `ggplot2` is perhaps more the way of the future?

Day 2 was mostly led by Ben Morris and Elliot Hauser. They went over the shell and version control, in particular git and github. We kept operating on our R work from Day 1. Ben also presented `make` and we continued to use our R work to demonstrate its power.

First task: tidy up the day 1 stuff. We create subdirectories: data, code, figures, results, prose. Then we move files into their respective homes. The second commit to this repo is a snapshot right after this. This was an exercise for the students, i.e. an opportunity for them to use many of the shell commands they just learned.

The file re-organization then requires us to visit the code files and prepend subdirectories, so the scripts still work. Making these incremental changes and documenting why they are necessary is how we started to demonstrate the power of version control. By this point Jenny had put up a public github repo (this one) and students were cloning / pulling from it. We got to see how `git pull` would not work if some of the student's local changes would be destroyed. We did not get into merging (for the students), so students were advised to discard their local changes, if they wanted to pull from this repo again.

Jenny and Ben continued to make changes to the repo, e.g. adding README.md and LICENSE, and providing various live demonstrations of using git and github to collaborate on this project. You will see some commits that are silly but are just us demonstrating something. Jenny and Ben purposely edited the same file and committed, so we could show them how merging worked.

After Ben had covered `make`, we set about making final preparations to some of our R scripts so they could be used in an automated pipeline managed by `make`. There are several commits related to these preparations. There is one remaining "gremlin" in which one of our automated reports is not successfully including a figure. Fixing this is on Jenny's to do list, it's probably merely a path issue relating to where the figures are being stored and sought.

Overall, we were quite pleased with how this all worked together. Day 1 stood well by itself and then that material gave us a great opportunity to apply everything we were learning on Day 2: file and directory manipulation from the command line, version control and collaboration, and automating a workflow with `make`.

PS: Jenny showed a couple slides from her UBC courses with helpful visuals for various R concepts. Students requested those and they are in `prose/slides.pdf`.
