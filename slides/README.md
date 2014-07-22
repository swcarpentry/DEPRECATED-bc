README
===

To create a new slideshow, copy `slides/pitch.html` and edit as necessary.

## Previewing the slides
From the top level directory for the `bc` repository, simply run:

`make site`

This will create the slides under `_sites/slides/`, where you can open
the relevant file in your browser to view it.

When you make an update to your slides and want to view them again, first
you will need to run

`make clean`

Before running

`make site`

again.

## Speaker Text
Text written inside the `<aside class="notes">` tag do not get rendered,
and instead serve as text for the speaker. This is important for other 
instructors who are less familiar with the material, and for those in 
other localities who will need to the text translated to their language.

## Assets
Images that appear in the slides should be stored within `imgs/slides/`

