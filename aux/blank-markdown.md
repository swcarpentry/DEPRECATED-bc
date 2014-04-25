---
layout: lesson
root: ../..
title: Title of Lesson
---
<div class="objectives" markdown="1">

#### Objectives
*   Describe the lesson's objectives in observable ways.
*   Not "understand X" (which isn't observable) but "Write a program that does X".
*   Aim for 3-6 points.
*   Keep the blank lines between the opening `<div ...>` and the closing `</div>` and these bullet points.

</div>

Write paragraphs of text here.
When you need to show input, output, and errors,
use `<div ...>` and triple tildes as shown below:

<div class="in" markdown="1">
~~~
$ this is the input
~~~
</div>
<div class="out" markdown="1">
~~~
this is the output
~~~
</div>
<div class="err" markdown="1">
~~~
error message
~~~
</div>

The `div`'s are needed because
Jekyll's Markdown processor will not let us put classes on code blocks.
We need these classes to be consistent with the HTML we produce for IPython Notebooks,
and because novices find examples much easier to read
when they can clearly distinguish input from output.

> #### Callout Boxes
> 
> Use a quoted block like this with a level-4 heading to write side notes.
> 
> These notes can span multiple paragraphs and include code blocks,
> but please try to keep them short.

Whenever you define a term,
include a [link](../../gloss.html#link)
to the `gloss.html` file in the root directory.
Note that this file doesn't exist in the repository:
it is produced during website compilation from `gloss.md`.

Images should be stored in the `img` directory below the lesson directory.
Please use SVG for diagrams,
since it scales better than raster formats like PNG or JPEG.
Please also include alternate text for accessibility aids and search engines:

~~~
<img src="img/filesystem.svg" alt="The Filesystem" />
~~~

<div class="keypoints" markdown="1">

#### Key Points
*   Every lesson should end with a summary of key points.
*   We will stitch these together to create reference guides for learners.
*   As with objectives, wrap this block in a `div` with the right style.

</div>

<div class="challenges" markdown="1">

#### Challenges

1.  Include a list of challenge exercises for learners at the end of the lesson.

2.  Put blank lines between the items so that the challenges will be spaced out in the final HTML document.

</div>
