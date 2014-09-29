---
layout: lesson
root: ../..
---

## Analyzing Patient Data


<<<<<<< HEAD
<div class="">
<p>We are studying inflammation in patients who have been given a new treatment for arthritis, and need to analyze the first dozen data sets. The data sets are stored in <a href="../../gloss.html#comma-separeted-values">comma-separated values</a> (CSV) format: each row holds information for a single patient, and the columns represent successive days. The first few rows of our first file look like this:</p>
<pre><code>0,0,1,3,1,2,4,7,8,3,3,3,10,5,7,4,7,7,12,18,6,13,11,11,7,7,4,6,8,8,4,4,5,7,3,4,2,3,0,0
0,1,2,1,2,1,3,2,2,6,10,11,5,9,4,4,7,16,8,6,18,4,12,5,12,7,11,5,11,3,3,5,4,4,5,5,1,1,0,1
0,1,1,3,3,2,6,2,5,9,5,7,4,5,4,15,5,11,9,10,19,14,12,17,7,12,11,7,4,2,10,5,4,2,2,3,2,2,1,1
0,0,2,0,4,2,2,1,6,7,10,7,9,13,8,8,15,10,10,7,17,4,4,7,6,15,6,4,9,11,3,5,6,3,3,4,2,3,2,1
0,1,1,3,3,1,3,5,2,4,4,7,6,5,3,10,8,10,6,17,9,14,9,7,13,9,12,6,7,7,9,6,3,2,2,4,2,0,1,1</code></pre>
</div>
=======
We are studying inflammation in patients who have been given a new treatment for arthritis,
and need to analyze the first dozen data sets.
The data sets are stored in [comma-separated values](../../gloss.html#comma-separeted-values) (CSV) format:
each row holds information for a single patient,
and the columns represent successive days.
The first few rows of our first file look like this:
>>>>>>> upstream/master

    0,0,1,3,1,2,4,7,8,3,3,3,10,5,7,4,7,7,12,18,6,13,11,11,7,7,4,6,8,8,4,4,5,7,3,4,2,3,0,0
    0,1,2,1,2,1,3,2,2,6,10,11,5,9,4,4,7,16,8,6,18,4,12,5,12,7,11,5,11,3,3,5,4,4,5,5,1,1,0,1
    0,1,1,3,3,2,6,2,5,9,5,7,4,5,4,15,5,11,9,10,19,14,12,17,7,12,11,7,4,2,10,5,4,2,2,3,2,2,1,1
    0,0,2,0,4,2,2,1,6,7,10,7,9,13,8,8,15,10,10,7,17,4,4,7,6,15,6,4,9,11,3,5,6,3,3,4,2,3,2,1
    0,1,1,3,3,1,3,5,2,4,4,7,6,5,3,10,8,10,6,17,9,14,9,7,13,9,12,6,7,7,9,6,3,2,2,4,2,0,1,1

<<<<<<< HEAD
<div class="">
<p>We want to:</p>
<ul>
<li>load that data into memory,</li>
<li>calculate the average inflammation per day across all patients, and</li>
<li>plot the result.</li>
</ul>
<p>To do all that, we'll have to learn a little bit about programming.</p>
</div>
=======
>>>>>>> upstream/master

We want to:

*   load that data into memory,
*   calculate the average inflammation per day across all patients, and
*   plot the result.

To do all that, we'll have to learn a little bit about programming.


<div class="objectives" markdown="1">
#### Objectives

*   Explain what a library is, and what libraries are used for.
*   Load a Python library and use the things it contains.
*   Read tabular data from a file into a program.
*   Assign values to variables.
*   Select individual values and subsections from data.
*   Perform operations on arrays of data.
*   Display simple graphs.
</div>

### Loading Data


<<<<<<< HEAD
<div class="">
<p>Words are useful, but what's more useful are the sentences and stories we use them to build. Similarly, while a lot of powerful tools are built into languages like Python, even more live in the <a href="../../gloss.html#library">libraries</a> they are used to build.</p>
<p>In order to load our inflammation data, we need to <a href="../../gloss.html#import">import</a> a library called NumPy that knows how to operate on matrices:</p>
</div>
=======
Words are useful,
but what's more useful are the sentences and stories we build with them.
Similarly,
while a lot of powerful tools are built into languages like Python,
even more live in the [libraries](../../gloss.html#library) they are used to build.

In order to load our inflammation data,
we need to [import](../../gloss.html#import) a library called NumPy. In general you should use this library if you want to do fancy things with numbers, especially if you have matrices. We can load NumPy using:
>>>>>>> upstream/master


<pre class="in"><code>import numpy</code></pre>


<<<<<<< HEAD
<div class="">
<p>Importing a library is like getting a piece of lab equipment out of a storage locker and setting it up on the bench. Once it's done, we can ask the library to read our data file for us:</p>
</div>
=======
Importing a library is like getting a piece of lab equipment out of a storage locker
and setting it up on the bench.
Once it's done,
we can ask the library to read our data file for us:
>>>>>>> upstream/master


<pre class="in"><code>numpy.loadtxt(fname=&#39;inflammation-01.csv&#39;, delimiter=&#39;,&#39;)</code></pre>

<div class="out"><pre class='out'><code>array([[ 0.,  0.,  1., ...,  3.,  0.,  0.],
       [ 0.,  1.,  2., ...,  1.,  0.,  1.],
       [ 0.,  1.,  1., ...,  2.,  1.,  1.],
       ..., 
       [ 0.,  1.,  1., ...,  1.,  1.,  1.],
       [ 0.,  0.,  0., ...,  0.,  2.,  0.],
       [ 0.,  0.,  1., ...,  1.,  1.,  0.]])</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>The expression <code>numpy.loadtxt(...)</code> is a <a href="../../gloss.html#function-call">function call</a> that asks Python to run the function <code>loadtxt</code> that belongs to the <code>numpy</code> library. This <a href="../../gloss.html#dotted-notation">dotted notation</a> is used everywhere in Python to refer to the parts of things as <code>whole.part</code>.</p>
<p><code>numpy.loadtxt</code> has two <a href="../../gloss.html#parameter">parameters</a>: the name of the file we want to read, and the <a href="../../gloss.html#delimiter">delimiter</a> that separates values on a line. These both need to be character strings (or <a href="../../gloss.html#string">strings</a> for short), so we put them in quotes.</p>
<p>When we are finished typing and press Shift+Enter, the notebook runs our command. Since we haven't told it to do anything else with the function's output, the notebook displays it. In this case, that output is the data we just loaded. By default, only a few rows and columns are shown (with <code>...</code> to omit elements when displaying big arrays). To save space, Python displays numbers as <code>1.</code> instead of <code>1.0</code> when there's nothing interesting after the decimal point.</p>
</div>


<div class="">
<p>Our call to <code>numpy.loadtxt</code> read our file, but didn't save the data in memory. To do that, we need to <a href="../../gloss.html#assignment">assign</a> the array to a <a href="../../gloss.html#variable">variable</a>. A variable is just a name for a value, such as <code>x</code>, <code>current_temperature</code>, or <code>subject_id</code>. We can create a new variable simply by assigning a value to it using <code>=</code>:</p>
</div>
=======
The expression `numpy.loadtxt(...)` is a [function call](../../gloss.html#function-call)
that asks Python to run the function `loadtxt` that belongs to the `numpy` library.
This [dotted notation](../../gloss.html#dotted-notation) is used everywhere in Python
to refer to the parts of things as `thing.component`.

`numpy.loadtxt` has two [parameters](../../gloss.html#parameter):
the name of the file we want to read,
and the [delimiter](../../gloss.html#delimiter) that separates values on a line.
These both need to be character strings (or [strings](../../gloss.html#string) for short),
so we put them in quotes.

When we are finished typing and press Shift+Enter,
the notebook runs our command.
Since we haven't told it to do anything else with the function's output,
the notebook displays it.
In this case,
that output is the data we just loaded.
By default,
only a few rows and columns are shown
(with `...` to omit elements when displaying big arrays).
To save space,
Python displays numbers as `1.` instead of `1.0`
when there's nothing interesting after the decimal point.


Our call to `numpy.loadtxt` read our file,
but didn't save the data in memory.
To do that,
we need to [assign](../../gloss.html#assignment) the array to a [variable](../../gloss.html#variable).
A variable is just a name for a value,
such as `x`, `current_temperature`, or `subject_id`. Python's variables must begin with a letter.
We can create a new variable simply by assigning a value to it using `=`:
>>>>>>> upstream/master


<pre class="in"><code>weight_kg = 55</code></pre>


Once a variable has a value, we can print it:


<pre class="in"><code>print weight_kg</code></pre>

<div class="out"><pre class='out'><code>55
</code></pre></div>


and do arithmetic with it:


<pre class="in"><code>print &#39;weight in pounds:&#39;, 2.2 * weight_kg</code></pre>

<div class="out"><pre class='out'><code>weight in pounds: 121.0
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>We can also change a variable's value by assigning it a new one:</p>
</div>
=======
We can also change a variable's value by assigning it a new one:
>>>>>>> upstream/master


<pre class="in"><code>weight_kg = 57.5
print &#39;weight in kilograms is now:&#39;, weight_kg</code></pre>

<div class="out"><pre class='out'><code>weight in kilograms is now: 57.5
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>As the example above shows, we can print several things at once by separating them with commas.</p>
<p>If we imagine the variable as a sticky note with a name written on it, assignment is like putting the sticky note on a particular value:</p>
</div>
=======
As the example above shows,
we can print several things at once by separating them with commas.
>>>>>>> upstream/master

If we imagine the variable as a sticky note with a name written on it,
assignment is like putting the sticky note on a particular value:


<img src="img/python-sticky-note-variables-01.svg" alt="Variables as Sticky Notes" />


<<<<<<< HEAD
<div class="">
<p>This means that assigning a value to one variable does <em>not</em> change the values of other variables. For example, let's store the subject's weight in pounds in a variable:</p>
</div>
=======
This means that assigning a value to one variable does *not* change the values of other variables.
For example,
let's store the subject's weight in pounds in a variable:
>>>>>>> upstream/master


<pre class="in"><code>weight_lb = 2.2 * weight_kg
print &#39;weight in kilograms:&#39;, weight_kg, &#39;and in pounds:&#39;, weight_lb</code></pre>

<div class="out"><pre class='out'><code>weight in kilograms: 57.5 and in pounds: 126.5
</code></pre></div>


<img src="img/python-sticky-note-variables-02.svg" alt="Creating Another Variable" />


and then change `weight_kg`:


<pre class="in"><code>weight_kg = 100.0
print &#39;weight in kilograms is now:&#39;, weight_kg, &#39;and weight in pounds is still:&#39;, weight_lb</code></pre>

<div class="out"><pre class='out'><code>weight in kilograms is now: 100.0 and weight in pounds is still: 126.5
</code></pre></div>


<img src="img/python-sticky-note-variables-03.svg" alt="Updating a Variable" />


<<<<<<< HEAD
<div class="">
<p>Since <code>weight_lb</code> doesn't &quot;remember&quot; where its value came from, it isn't automatically updated when <code>weight_kg</code> changes. This is different from the way spreadsheets work.</p>
<p>Now that we know how to assign things to variables, let's re-run <code>numpy.loadtxt</code> and save its result:</p>
</div>
=======
Since `weight_lb` doesn't "remember" where its value came from,
it isn't automatically updated when `weight_kg` changes.
This is different from the way spreadsheets work.

Now that we know how to assign things to variables,
let's re-run `numpy.loadtxt` and save its result:
>>>>>>> upstream/master


<pre class="in"><code>data = numpy.loadtxt(fname=&#39;inflammation-01.csv&#39;, delimiter=&#39;,&#39;)</code></pre>


<<<<<<< HEAD
<div class="">
<p>This statement doesn't produce any output because assignment doesn't display anything. If we want to check that our data has been loaded, we can print the variable's value:</p>
</div>
=======
This statement doesn't produce any output because assignment doesn't display anything.
If we want to check that our data has been loaded,
we can print the variable's value:
>>>>>>> upstream/master


<pre class="in"><code>print data</code></pre>

<div class="out"><pre class='out'><code>[[ 0.  0.  1. ...,  3.  0.  0.]
 [ 0.  1.  2. ...,  1.  0.  1.]
 [ 0.  1.  1. ...,  2.  1.  1.]
 ..., 
 [ 0.  1.  1. ...,  1.  1.  1.]
 [ 0.  0.  0. ...,  0.  2.  0.]
 [ 0.  0.  1. ...,  1.  1.  0.]]
</code></pre></div>


<<<<<<< HEAD
<div class="challenges">
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Draw diagrams showing what variables refer to what values after each statement in the following program:</p>
<pre class="sourceCode python"><code class="sourceCode python">mass = <span class="fl">47.5</span>
age = <span class="dv">122</span>
mass = mass * <span class="fl">2.0</span>
age = age - <span class="dv">20</span></code></pre></li>
<li><p>Again, draw diagrams showing what variables refer to what values after each statement. What does the following program print out?</p>
<pre class="sourceCode python"><code class="sourceCode python">item1, value1 = <span class="st">&#39;Disease&#39;</span>, <span class="st">&#39;Arthritis&#39;</span>
item2, value2 = value1, item1
<span class="dt">print</span> item2, value2</code></pre></li>
</ol>
=======
<div class="challenges" markdown="1">
#### Challenges

1.  Draw diagrams showing what variables refer to what values after each statement in the following program:

    ~~~python
    mass = 47.5
    age = 122
    mass = mass * 2.0
    age = age - 20
    ~~~

1.  What does the following program print out?

    ~~~python
    first, second = 'Grace', 'Hopper'
    third, fourth = second, first
    print third, fourth
    ~~~
>>>>>>> upstream/master
</div>

### Manipulating Data


<<<<<<< HEAD
<div class="">
<p>Now that our data is in memory, we can start doing things with it. First, let's ask what <a href="../../gloss.html#data-type">type</a> of thing <code>data</code> refers to:</p>
</div>
=======
Now that our data is in memory,
we can start doing things with it.
First,
let's ask what [type](../../gloss.html#data-type) of thing `data` refers to:
>>>>>>> upstream/master


<pre class="in"><code>print type(data)</code></pre>

<div class="out"><pre class='out'><code>&lt;type &#39;numpy.ndarray&#39;&gt;
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>The output tells us that <code>data</code> currently refers to an N-dimensional array created by the NumPy library. We can see what its <a href="../../gloss.html#shape">shape</a> is like this:</p>
</div>
=======
The output tells us that `data` currently refers to an N-dimensional array created by the NumPy library.
We can see what its [shape](../../gloss.html#shape) is like this:
>>>>>>> upstream/master


<pre class="in"><code>print data.shape</code></pre>

<div class="out"><pre class='out'><code>(60, 40)
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>This tells us that <code>data</code> has 60 rows and 40 columns. <code>data.shape</code> is a <a href="../../gloss.html#member">member</a> of <code>data</code>, i.e., a value that is stored as part of a larger value. We use the same dotted notation for the members of values that we use for the functions in libraries because they have the same part-and-whole relationship.</p>
</div>


<div class="">
<p>If we want to get a single value from the matrix, we must provide an <a href="../../gloss.html#index">index</a> in square brackets, just as we do in math:</p>
</div>
=======
This tells us that `data` has 60 rows and 40 columns.
`data.shape` is a [member](../../gloss.html#member) of `data`,
i.e.,
a value that is stored as part of a larger value.
We use the same dotted notation for the members of values
that we use for the functions in libraries
because they have the same part-and-whole relationship.


If we want to get a single value from the matrix,
we must provide an [index](../../gloss.html#index) in square brackets,
just as we do in math:
>>>>>>> upstream/master


<pre class="in"><code>print &#39;first value in data:&#39;, data[0, 0]</code></pre>

<div class="out"><pre class='out'><code>first value in data: 0.0
</code></pre></div>


<pre class="in"><code>print &#39;middle value in data:&#39;, data[30, 20]</code></pre>

<div class="out"><pre class='out'><code>middle value in data: 13.0
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>The expression <code>data[30, 20]</code> may not surprise you, but <code>data[0, 0]</code> might. Programming languages like Fortran and MATLAB start counting at 1, because that's what human beings have done for thousands of years. Languages in the C family (including C++, Java, Perl, and Python) count from 0 because that's simpler for computers to do. As a result, if we have an M×N array in Python, its indices go from 0 to M-1 on the first axis and 0 to N-1 on the second. It takes a bit of getting used to, but one way to remember the rule is that the index is how many steps we have to take from the start to get the item we want.</p>
<blockquote>
<h4>In the Corner</h4>
<p>What may also surprise you is that when Python displays an array, it shows the element with index <code>[0, 0]</code> in the upper left corner rather than the lower left. This is consistent with the way mathematicians draw matrices, but different from the Cartesian coordinates. The indices are (row, column) instead of (column, row) for the same reason.</p>
</blockquote>
</div>
=======
The expression `data[30, 20]` may not surprise you,
but `data[0, 0]` might.
Programming languages like Fortran and MATLAB start counting at 1,
because that's what human beings have done for thousands of years.
Languages in the C family (including C++, Java, Perl, and Python) count from 0
because that's simpler for computers to do.
As a result,
if we have an M&times;N array in Python,
its indices go from 0 to M-1 on the first axis
and 0 to N-1 on the second.
It takes a bit of getting used to,
but one way to remember the rule is that
the index is how many steps we have to take from the start to get the item we want.
>>>>>>> upstream/master

> #### In the Corner
>
> What may also surprise you is that when Python displays an array,
> it shows the element with index `[0, 0]` in the upper left corner
> rather than the lower left.
> This is consistent with the way mathematicians draw matrices,
> but different from the Cartesian coordinates.
> The indices are (row, column) instead of (column, row) for the same reason, 
> which can be confusing when plotting data.

<<<<<<< HEAD
<div class="">
<p>An index like <code>[30, 20]</code> selects a single element of an array, but we can select whole sections as well. For example, we can select the first ten days (columns) of values for the first four (rows) patients like this:</p>
</div>
=======

An index like `[30, 20]` selects a single element of an array,
but we can select whole sections as well.
For example,
we can select the first ten days (columns) of values
for the first four (rows) patients like this:
>>>>>>> upstream/master


<pre class="in"><code>print data[0:4, 0:10]</code></pre>

<div class="out"><pre class='out'><code>[[ 0.  0.  1.  3.  1.  2.  4.  7.  8.  3.]
 [ 0.  1.  2.  1.  2.  1.  3.  2.  2.  6.]
 [ 0.  1.  1.  3.  3.  2.  6.  2.  5.  9.]
 [ 0.  0.  2.  0.  4.  2.  2.  1.  6.  7.]]
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>The <a href="../../gloss.html#slice">slice</a> <code>0:4</code> means, &quot;Start at index 0 and go up to, but not including, index 4.&quot; Again, the up-to-but-not-including takes a bit of getting used to, but the rule is that the difference between the upper and lower bounds is the number of values in the slice.</p>
<p>We don't have to start slices at 0:</p>
</div>
=======
The [slice](../../gloss.html#slice) `0:4` means,
"Start at index 0 and go up to, but not including, index 4."
Again,
the up-to-but-not-including takes a bit of getting used to,
but the rule is that the difference between the upper and lower bounds is the number of values in the slice.

We don't have to start slices at 0:
>>>>>>> upstream/master


<pre class="in"><code>print data[5:10, 0:10]</code></pre>

<div class="out"><pre class='out'><code>[[ 0.  0.  1.  2.  2.  4.  2.  1.  6.  4.]
 [ 0.  0.  2.  2.  4.  2.  2.  5.  5.  8.]
 [ 0.  0.  1.  2.  3.  1.  2.  3.  5.  3.]
 [ 0.  0.  0.  3.  1.  5.  6.  5.  5.  8.]
 [ 0.  1.  1.  2.  1.  3.  5.  3.  5.  8.]]
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>and we don't have to take all the values in the slice---if we provide a <a href="../../gloss.html#stride">stride</a>, Python takes values spaced that far apart:</p>
</div>


<pre class="in"><code>print data[0:10:3, 0:10:2]</code></pre>

<div class="out"><pre class='out'><code>[[ 0.  1.  1.  4.  8.]
 [ 0.  2.  4.  2.  6.]
 [ 0.  2.  4.  2.  5.]
 [ 0.  1.  1.  5.  5.]]
</code></pre></div>


<div class="">
<p>Here, we have taken rows 0, 3, 6, and 9, and columns 0, 2, 4, 6, and 8. (Again, we always include the lower bound, but stop when we reach or cross the upper bound.)</p>
</div>


<div class="">
<p>We also don't have to include the upper and lower bound on the slice. If we don't include the lower bound, Python uses 0 by default; if we don't include the upper, the slice runs to the end of the axis, and if we don't include either (i.e., if we just use ':' on its own), the slice includes everything:</p>
</div>
=======
We also don't have to include the upper and lower bound on the slice.
If we don't include the lower bound,
Python uses 0 by default;
if we don't include the upper,
the slice runs to the end of the axis,
and if we don't include either
(i.e., if we just use ':' on its own),
the slice includes everything:
>>>>>>> upstream/master


<pre class="in"><code>small = data[:3, 36:]
print &#39;small is:&#39;
print small</code></pre>

<div class="out"><pre class='out'><code>small is:
[[ 2.  3.  0.  0.]
 [ 1.  1.  0.  1.]
 [ 2.  2.  1.  1.]]
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>Arrays also know how to perform common mathematical operations on their values. If we want to find the average inflammation for all patients on all days, for example, we can just ask the array for its mean value</p>
</div>
=======
Arrays also know how to perform common mathematical operations on their values.
If we want to find the average inflammation for all patients on all days,
for example,
we can just ask the array for its mean value
>>>>>>> upstream/master


<pre class="in"><code>print data.mean()</code></pre>

<div class="out"><pre class='out'><code>6.14875
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p><code>mean</code> is a <a href="../../gloss.html#method">method</a> of the array, i.e., a function that belongs to it in the same way that the member <code>shape</code> does. If variables are nouns, methods are verbs: they are what the thing in question knows how to do. This is why <code>data.shape</code> doesn't need to be called (it's just a thing) but <code>data.mean()</code> does (it's an action). It is also why we need empty parentheses for <code>data.mean()</code>: even when we're not passing in any parameters, parentheses are how we tell Python to go and do something for us.</p>
<p>NumPy arrays have lots of useful methods:</p>
</div>
=======
`mean` is a [method](../../gloss.html#method) of the array,
i.e.,
a function that belongs to it
in the same way that the member `shape` does.
If variables are nouns, methods are verbs:
they are what the thing in question knows how to do.
This is why `data.shape` doesn't need to be called
(it's just a thing)
but `data.mean()` does
(it's an action).
It is also why we need empty parentheses for `data.mean()`:
even when we're not passing in any parameters,
parentheses are how we tell Python to go and do something for us.

NumPy arrays have lots of useful methods:
>>>>>>> upstream/master


<pre class="in"><code>print &#39;maximum inflammation:&#39;, data.max()
print &#39;minimum inflammation:&#39;, data.min()
print &#39;standard deviation:&#39;, data.std()</code></pre>

<div class="out"><pre class='out'><code>maximum inflammation: 20.0
minimum inflammation: 0.0
standard deviation: 4.61383319712
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>When analyzing data, though, we often want to look at partial statistics, such as the maximum value per patient or the average value per day. One way to do this is to select the data we want to create a new temporary array, then ask it to do the calculation:</p>
</div>
=======
When analyzing data,
though,
we often want to look at partial statistics,
such as the maximum value per patient
or the average value per day.
One way to do this is to select the data we want to create a new temporary array,
then ask it to do the calculation:
>>>>>>> upstream/master


<pre class="in"><code>patient_0 = data[0, :] # 0 on the first axis, everything on the second
print &#39;maximum inflammation for patient 0:&#39;, patient_0.max()</code></pre>

<div class="out"><pre class='out'><code>maximum inflammation for patient 0: 18.0
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>We don't actually need to store the row in a variable of its own. Instead, we can combine the selection and the method call:</p>
</div>
=======
We don't actually need to store the row in a variable of its own.
Instead, we can combine the selection and the method call:
>>>>>>> upstream/master


<pre class="in"><code>print &#39;maximum inflammation for patient 2:&#39;, data[2, :].max()</code></pre>

<div class="out"><pre class='out'><code>maximum inflammation for patient 2: 19.0
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>What if we need the maximum inflammation for <em>all</em> patients, or the average for each day? As the diagram below shows, we want to perform the operation across an axis:</p>
</div>
=======
What if we need the maximum inflammation for *all* patients,
or the average for each day?
As the diagram below shows,
we want to perform the operation across an axis:
>>>>>>> upstream/master


<img src="img/python-operations-across-axes.svg" alt="Operations Across Axes" />


<<<<<<< HEAD
<div class="">
<p>To support this, most array methods allow us to specify the axis we want to work on. If we ask for the average across axis 0, we get:</p>
</div>
=======
To support this,
most array methods allow us to specify the axis we want to work on.
If we ask for the average across axis 0,
we get:
>>>>>>> upstream/master


<pre class="in"><code>print data.mean(axis=0)</code></pre>

<div class="out"><pre class='out'><code>[  0.           0.45         1.11666667   1.75         2.43333333   3.15
   3.8          3.88333333   5.23333333   5.51666667   5.95         5.9
   8.35         7.73333333   8.36666667   9.5          9.58333333
  10.63333333  11.56666667  12.35        13.25        11.96666667
  11.03333333  10.16666667  10.           8.66666667   9.15         7.25
   7.33333333   6.58333333   6.06666667   5.95         5.11666667   3.6
   3.3          3.56666667   2.48333333   1.5          1.13333333
   0.56666667]
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>As a quick check, we can ask this array what its shape is:</p>
</div>
=======
As a quick check,
we can ask this array what its shape is:
>>>>>>> upstream/master


<pre class="in"><code>print data.mean(axis=0).shape</code></pre>

<div class="out"><pre class='out'><code>(40,)
</code></pre></div>


<<<<<<< HEAD
<div class="">
<p>The expression <code>(40,)</code> tells us we have an N×1 vector, so this is the average inflammation per day for all patients. If we average across axis 1, we get:</p>
</div>
=======
The expression `(40,)` tells us we have an N&times;1 vector,
so this is the average inflammation per day for all patients.
If we average across axis 1, we get:
>>>>>>> upstream/master


<pre class="in"><code>print data.mean(axis=1)</code></pre>

<div class="out"><pre class='out'><code>[ 5.45   5.425  6.1    5.9    5.55   6.225  5.975  6.65   6.625  6.525
  6.775  5.8    6.225  5.75   5.225  6.3    6.55   5.7    5.85   6.55
  5.775  5.825  6.175  6.1    5.8    6.425  6.05   6.025  6.175  6.55
  6.175  6.35   6.725  6.125  7.075  5.725  5.925  6.15   6.075  5.75
  5.975  5.725  6.3    5.9    6.75   5.925  7.225  6.15   5.95   6.275  5.7
  6.1    6.825  5.975  6.725  5.7    6.25   6.4    7.05   5.9  ]
</code></pre></div>


which is the average inflammation per patient across all days.


<div class="challenges" markdown="1">
#### Challenges

<<<<<<< HEAD
<div class="challenges">
<h4 id="challenges">Challenges</h4>
<p>A subsection of an array is called a <a href="../../gloss.html#slice">slice</a>. We can take slices of character strings as well:</p>
=======
A subsection of an array is called a [slice](../../gloss.html#slice).
We can take slices of character strings as well:
>>>>>>> upstream/master
</div>


<pre class="in"><code>element = &#39;oxygen&#39;
print &#39;first three characters:&#39;, element[0:3]
print &#39;last three characters:&#39;, element[3:6]</code></pre>

<div class="out"><pre class='out'><code>first three characters: oxy
last three characters: gen
</code></pre></div>


<<<<<<< HEAD
<div class="challenges">
<ol style="list-style-type: decimal">
<li><p>What is the value of <code>element[:4]</code>? What about <code>element[4:]</code>? Or <code>element[:]</code>?</p></li>
<li><p>What is <code>element[-1]</code>? What is <code>element[-2]</code>? Given those answers, explain what <code>element[1:-1]</code> does.</p></li>
<li><p>The expression <code>element[3:3]</code> produces an <a href="../../gloss.html#empty-string">empty string</a>, i.e., a string that contains no characters. If <code>data</code> holds our array of patient data, what does <code>data[3:3, 4:4]</code> produce? What about <code>data[3:3, :]</code>?</p></li>
</ol>
=======
<div class="challenges" markdown="1">
1.  What is the value of `element[:4]`?
    What about `element[4:]`?
    Or `element[:]`?

1.  What is `element[-1]`?
    What is `element[-2]`?
    Given those answers,
    explain what `element[1:-1]` does.

1.  The expression `element[3:3]` produces an [empty string](../../gloss.html#empty-string),
    i.e., a string that contains no characters.
    If `data` holds our array of patient data,
    what does `data[3:3, 4:4]` produce?
    What about `data[3:3, :]`?
>>>>>>> upstream/master
</div>

### Plotting


<<<<<<< HEAD
<div class="">
<p>The mathematician Richard Hamming once said, &quot;The purpose of computing is insight, not numbers,&quot; and the best way to develop insight is often to visualize data. Visualization deserves an entire lecture (or course) of its own, but we can explore a few features of Python's <code>matplotlib</code> here. First, let's tell the IPython Notebook that we want our plots displayed inline, rather than in a separate viewing window:</p>
</div>
=======
The mathematician Richard Hamming once said,
"The purpose of computing is insight, not numbers,"
and the best way to develop insight is often to visualize data.
Visualization deserves an entire lecture (or course) of its own,
but we can explore a few features of Python's `matplotlib` here. While there is no "official" plotting library, this package is the de facto standard.
First,
let's tell the IPython Notebook that we want our plots displayed inline,
rather than in a separate viewing window:
>>>>>>> upstream/master


<pre class="in"><code>%matplotlib inline</code></pre>


<<<<<<< HEAD
<div class="">
<p>The <code>%</code> at the start of the line signals that this is a command for the notebook, rather than a statement in Python. Next, we will import the <code>pyplot</code> module from <code>matplotlib</code> and use two of its functions to create and display a heat map of our data:</p>
</div>
=======
The `%` at the start of the line signals that this is a command for the notebook,
rather than a statement in Python.
Next,
we will import the `pyplot` module from `matplotlib`
and use two of its functions to create and display a heat map of our data:
>>>>>>> upstream/master


<pre class="in"><code>from matplotlib import pyplot
pyplot.imshow(data)
pyplot.show()</code></pre>

<div class="out">
<img src="../../novice/python/01-numpy_files/novice/python/01-numpy_71_0.png">
</div>


<<<<<<< HEAD
<div class="">
<p>Blue regions in this heat map are low values, while red shows high values. As we can see, inflammation rises and falls over a 40-day period. Let's take a look at the average inflammation over time:</p>
</div>
=======
Blue regions in this heat map are low values, while red shows high values.
As we can see,
inflammation rises and falls over a 40-day period.
Let's take a look at the average inflammation over time:
>>>>>>> upstream/master


<pre class="in"><code>ave_inflammation = data.mean(axis=0)
pyplot.plot(ave_inflammation)
pyplot.show()</code></pre>

<div class="out">
<img src="../../novice/python/01-numpy_files/novice/python/01-numpy_73_0.png">
</div>


<<<<<<< HEAD
<div class="">
<p>Here, we have put the average per day across all patients in the variable <code>ave_inflammation</code>, then asked <code>pyplot</code> to create and display a line graph of those values. The result is roughly a linear rise and fall, which is suspicious: based on other studies, we expect a sharper rise and slower fall. Let's have a look at two other statistics:</p>
</div>
=======
Here,
we have put the average per day across all patients in the variable `ave_inflammation`,
then asked `pyplot` to create and display a line graph of those values.
The result is roughly a linear rise and fall,
which is suspicious:
based on other studies,
we expect a sharper rise and slower fall.
Let's have a look at two other statistics:
>>>>>>> upstream/master


<pre class="in"><code>print &#39;maximum inflammation per day&#39;
pyplot.plot(data.max(axis=0))
pyplot.show()

print &#39;minimum inflammation per day&#39;
pyplot.plot(data.min(axis=0))
pyplot.show()</code></pre>

<div class="out"><pre class='out'><code>maximum inflammation per day
</code></pre>
<img src="../../novice/python/01-numpy_files/novice/python/01-numpy_75_1.png">
<pre class='out'><code>minimum inflammation per day
</code></pre>
<img src="../../novice/python/01-numpy_files/novice/python/01-numpy_75_3.png">
</div>


<<<<<<< HEAD
<div class="">
<p>The maximum value rises and falls perfectly smoothly, while the minimum seems to be a step function. Neither result seems particularly likely, so either there's a mistake in our calculations or something is wrong with our data.</p>
</div>
=======
The maximum value rises and falls perfectly smoothly,
while the minimum seems to be a step function.
Neither result seems particularly likely,
so either there's a mistake in our calculations
or something is wrong with our data.

>>>>>>> upstream/master

<div class="challenges" markdown="1">
#### Challenges

<<<<<<< HEAD
<div class="challenges">
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Why do all of our plots stop just short of the upper end of our graph? Why are the vertical lines in our plot of the minimum inflammation per day not vertical?</p></li>
<li><p>Create a plot showing the standard deviation of the inflammation data for each day across all patients.</p></li>
</ol>
=======
1.  Why do all of our plots stop just short of the upper end of our graph?
    Why are the vertical lines in our plot of the minimum inflammation per day not perfectly vertical?

2.  Create a plot showing the standard deviation of the inflammation data for each day across all patients.
>>>>>>> upstream/master
</div>

### Wrapping Up


<<<<<<< HEAD
<div class="">
<p>It's very common to create an <a href="../../gloss.html#alias">alias</a> for a library when importing it in order to reduce the amount of typing we have to do. Here are our three plots side by side using aliases for <code>numpy</code> and <code>pyplot</code>:</p>
</div>
=======
It's very common to create an [alias](../../gloss.html#alias) for a library when importing it
in order to reduce the amount of typing we have to do.
Here are our three plots side by side using aliases for `numpy` and `pyplot`:
>>>>>>> upstream/master


<pre class="in"><code>import numpy as np
from matplotlib import pyplot as plt

data = np.loadtxt(fname=&#39;inflammation-01.csv&#39;, delimiter=&#39;,&#39;)

plt.figure(figsize=(10.0, 3.0))

plt.subplot(1, 3, 1)
plt.ylabel(&#39;average&#39;)
plt.plot(data.mean(0))

plt.subplot(1, 3, 2)
plt.ylabel(&#39;max&#39;)
plt.plot(data.max(0))

plt.subplot(1, 3, 3)
plt.ylabel(&#39;min&#39;)
plt.plot(data.min(0))

plt.tight_layout()
plt.show()</code></pre>

<div class="out">
<img src="../../novice/python/01-numpy_files/novice/python/01-numpy_80_0.png">
</div>


<<<<<<< HEAD
<div class="">
<p>The first two lines re-load our libraries as <code>np</code> and <code>plt</code>, which are the aliases most Python programmers use. The call to <code>loadtxt</code> reads our data, and the rest of the program tells the plotting library how large we want the figure to be, that we're creating three sub-plots, what to draw for each one, and that we want a tight layout. (Perversely, if we leave out that call to <code>plt.tight_layout()</code>, the graphs will actually be squeezed together more closely.)</p>
</div>
=======
The first two lines re-load our libraries as `np` and `plt`,
which are the aliases most Python programmers use.
The call to `loadtxt` reads our data,
and the rest of the program tells the plotting library
how large we want the figure to be,
that we're creating three sub-plots,
what to draw for each one,
and that we want a tight layout.
(Perversely,
if we leave out that call to `plt.tight_layout()`,
the graphs will actually be squeezed together more closely.)

>>>>>>> upstream/master

<div class="challenges" markdown="1">
#### Challenges

<<<<<<< HEAD
<div class="challenges">
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li>Modify the program to display the three plots on top of one another instead of side by side.</li>
</ol>
=======
1.  Modify the program to display the three plots on top of one another instead of side by side.
>>>>>>> upstream/master
</div>


<div class="keypoints" markdown="1">
#### Key Points

*   Import a library into a program using `import libraryname`.
*   Use the `numpy` library to work with arrays in Python.
*   Use `variable = value` to assign a value to a variable in order to record it in memory.
*   Variables are created on demand whenever a value is assigned to them.
*   Use `print something` to display the value of `something`.
*   The expression `array.shape` gives the shape of an array.
*   Use `array[x, y]` to select a single element from an array.
*   Array indices start at 0, not 1.
*   Use `low:high` to specify a slice that includes the indices from `low` to `high-1`.
*   All the indexing and slicing that works on arrays also works on strings.
*   Use `# some kind of explanation` to add comments to programs.
*   Use `array.mean()`, `array.max()`, and `array.min()` to calculate simple statistics.
*   Use `array.mean(axis=0)` or `array.mean(axis=1)` to calculate statistics across the specified axis.
*   Use the `pyplot` library from `matplotlib` for creating simple visualizations.
</div>


<<<<<<< HEAD
<div class="">
<h4 id="next-steps">Next Steps</h4>
<p>Our work so far has convinced us that something's wrong with our first data file. We would like to check the other 11 the same way, but typing in the same commands repeatedly is tedious and error-prone. Since computers don't get bored (that we know of), we should create a way to do a complete analysis with a single command, and then figure out how to repeat that step once for each file. These operations are the subjects of the next two lessons.</p>
</div>
=======
#### Next Steps

Our work so far has convinced us that something's wrong with our first data file.
We would like to check the other 11 the same way,
but typing in the same commands repeatedly is tedious and error-prone.
Since computers don't get bored (that we know of),
we should create a way to do a complete analysis with a single command,
and then figure out how to repeat that step once for each file.
These operations are the subjects of the next two lessons.
>>>>>>> upstream/master
