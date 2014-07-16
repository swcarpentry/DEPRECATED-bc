---
layout: lesson
root: ../..
---

## Analyzing Patient Data


We are studying inflammation in patients who have been given a new treatment for arthritis,
and need to analyze the first dozen data sets.
The data sets are stored in [CSV](../../gloss.html#comma-separated-values) format:
each row holds information for a single patient,
and the columns represent successive days.
The first few rows of our first file, 
`inflammation-01.csv`, look like this:

~~~
0,0,1,3,1,2,4,7,8,3,3,3,10,5,7,4,7,7,12,18,6,13,11,11,7,7,4,6,8,8,4,4,5,7,3,4,2,3,0,0
0,1,2,1,2,1,3,2,2,6,10,11,5,9,4,4,7,16,8,6,18,4,12,5,12,7,11,5,11,3,3,5,4,4,5,5,1,1,0,1
0,1,1,3,3,2,6,2,5,9,5,7,4,5,4,15,5,11,9,10,19,14,12,17,7,12,11,7,4,2,10,5,4,2,2,3,2,2,1,1
0,0,2,0,4,2,2,1,6,7,10,7,9,13,8,8,15,10,10,7,17,4,4,7,6,15,6,4,9,11,3,5,6,3,3,4,2,3,2,1
0,1,1,3,3,1,3,5,2,4,4,7,6,5,3,10,8,10,6,17,9,14,9,7,13,9,12,6,7,7,9,6,3,2,2,4,2,0,1,1
~~~
{:class="out"}



We want to:

* load that data into memory,
* calculate the average inflammation per day across all patients, and
* plot the result.

To do all that, we'll have to learn a little bit about programming.

<div class="objectives" markdown="1">
#### Objectives
* Learn about Matlab arrays.
* Read tabular data from a file into a program.
* Assign values to variables.
* Select individual values and subsections from data.
* Perform operations on arrays of data.
* Display simple graphs.
</div>

### Loading Data

Reading data from files and writing data to them
are essential tasks in scientific computing, 
and admittedly, something that we'd rather not spend a lot of time 
thinking about. Fortunately, Matlab comes with a number of high-level
tools to do these things 
efficiently, sparing us the grisly detail.

To fetch the data from our CSV file into Matlab, type following 
command into the Matlab shell, and press Enter:

~~~
csvread('inflammation-01.csv')
~~~
{:class="in"}

You should see a wall of numbers on the screen---these are the values
from the CSV file.
It can sometimes
be useful to see the output from Matlab commands, but it is often not.
To suppress the
output, simply put a semicolon at the end of your command:

~~~
csvread('inflammation-01.csv');
~~~
{:class="in"}

The expression `csvread(...)` is a 
[function call](../../gloss.html#function-call). 
Functions generally need [parameters](../../gloss.html#parameter)
to run.
In the case of the `csvread` function, we need to provide a single
parameter: the name of the file we want to read data from. This 
parameter needs to be a character string or 
[string](../../gloss.html#string), so we put it in quotes.

Our call to `csvread` read our file, and printed the data inside 
to the screen. But we still have no way to modify those values
or compute with them. To do that, we need to assign the array to a
[variable](../../gloss.html#variable).

~~~
patient_data = csvread('inflammation-01.csv');
~~~
{:class="in"}

A variable is just a name for a piece of data or *value*. 
Variable names must begin with a letter, and can contain 
numbers or underscores. Examples of valid variable names are 
`x`, `f_0` or `current_temperature`.

We can create a new variable simply by assigning a value to it using
`=`:

~~~
weight_kg = 55;
~~~
{:class="in"}

Once a variable has a value, we can print it using the `disp` function:

~~~
disp(weight_kg);
~~~
{:class="in"}

and do arithmetic with it:

~~~
weight_in_pounds = 2.2 * weight_kg;
disp(['Weight in pounds: ', num2str(weight_in_pounds)]);
~~~
{:class="in"}

~~~
Weight in pounds: 121.0
~~~
{:class="out"}

The `disp` function takes a single parameter -- the value to print. To 
print more than one value on a single line, we could print an *array*
of values. All values in this array need to be the same type. So, if 
we want to print a string and a numerical value together, we *have* to
convert that numerical value to a string with the `num2str` function.


#### Understanding Assignment

If we imagine the variable as a sticky note with a name written on 
it, assignment is like putting the sticky note on a particular value: 


<div>
  <img src="img/matlab-sticky-note-variables-01.svg" alt="Variables as Sticky Notes" />
</div>



Assigning a value to one variable does not change the values of other 
variables.
For example,

~~~
weight_kg = 57.5;
weight_lb = 2.2 * weight_in__kg;
disp(['Weight in kg: ', num2str(weight_kg); 'Weight in pounds: ', num2str(weight_lb)]);
~~~
{:class="in"}

~~~
Weight in kg: 57.5
Weight in pounds: 126.5
~~~
{:class="out"}


<div>
  <img src="img/matlab-sticky-note-variables-02.svg" alt="Creating another variable" />
</div>



Let's update the value of one of our variable, and print the values
of both:

~~~
weight_kg = 100;
disp(['Weight in kg: ', num2str(weight_kg); 'Weight in pounds: ', 
num2str(weight_lb)]);
~~~
{:class="in"}

~~~
Weight in kg: 100
Weight in pounds: 126.5
~~~
{:class="out"}


<div>
  <img src="img/matlab-sticky-note-variables-03.svg" alt="Updating one variable" />
</div>

Since `weight_lb` doesn't "remember" where its value came from, it isnt
automatically updated when `weight_kg` changes. This is important to
remember, and different from the way spreadsheets work.

Now that we know how to assign things to variables, let's re-run
`csvread` and save its result.

~~~
patient_data = csvread('inflammation-01.csv');
~~~
{:class="in"}

Matlab provides a command
to list all variables that have been assigned data:

~~~
who
~~~
{:class="in"}

~~~
Variables in the current scope:


patient_data
~~~
{:class="out"}


### Manipulating Data

Now that our data is in memory, we can start doing things with it.
First, let's find out its size or [shape](../../gloss.html#shape):

~~~
size(patient_data)
~~~
{:class="in"}


~~~
ans =
    
    60 40
~~~
{:class="out"}

The output tells us that the variable `patient_data` 
refers to a table of values
that has 60 rows and 40 columns.

Matlab stores *all* data in the form of arrays. For example:

* Numbers, or *scalars* are arrays of zero dimensions, as are single 
characters,
* Lists of numbers, or *vectors* are arrays of one dimension,
* Tables of numbers, or *matrices* are arrays of two dimensions,
* Even character strings, like sentences, are stored as an "array
of characters".

Normally, Matlab arrays can't store elements of different data types. For
instance, a Matlab array can't store both an `int` and a `char`. To do that, 
you have to use a [Cell Array](http://www.mathworks.com/help/matlab/cell-arrays.html). 

We can use the `class` function to find out what kind of data lives
inside an array:

~~~
class(patient_data)
~~~
{:class="in"}

~~~
ans = double
~~~
{:class="out"}

This output tells us that `patient_data` refers to an array of 
double precision floating-point numbers.

If we want to get a single value from the matrix, we must provide
an [index](../../gloss.html#index) in brackets, just as we do in math:

~~~
patient_data(1, 1)
~~~
{:class="in"}

~~~
ans = 0
~~~
{:class="out"}

This means that the value sitting at the `(1, 1)` position of the table,
i.e., the first row and first column, is `0`,

~~~
patient_data(30, 20)
~~~
{:class="in"}

~~~
ans = 16
~~~
{:class="out"}

and the value corresponding to the 30th row and 20th column, is `16`.


What may surprise you is that when Matlab displays an array, it shows
the element with index `(1, 1)` position in the upper left corner rather
than the lower left. This is consistent with the way mathematicians draw 
matrices, but different from the Cartesian coordinates.
The indices are (row, column) instead of (column, row) for the same
reason.

An index like `(30, 20)` selects a single element of 
an array, but we can select whole sections as well. For example,
we can select the first ten days (columns) of values for the first
four (rows) patients like this:

~~~
patient_data(1:4, 1:10)
~~~
{:class="in"}

~~~
ans =

   0   0   1   3   1   2   4   7   8   3
   0   1   2   1   2   1   3   2   2   6
   0   1   1   3   3   2   6   2   5   9
   0   0   2   0   4   2   2   1   6   7
~~~
{:class="out"}

The [slice](../../gloss.html#slice) `(1:4)` means, "Start at index
1 and go up to 4". 

We don't have to start slices at 1:

~~~
patient_data(6:10, 1:10)
~~~
{:class="in"}

~~~
ans =

   0   0   1   2   2   4   2   1   6   4
   0   0   2   2   4   2   2   5   5   8
   0   0   1   2   3   1   2   3   5   3
   0   0   0   3   1   5   6   5   5   8
   0   1   1   2   1   3   5   3   5   8
~~~
{:class="out"}

and we don't have to take all the values in the slice---if we provide
a [stride](../../gloss.html#stride),

~~~
patient_data(1:3:10, 1:2:10)
~~~
{:class="in"}

~~~
ans =

   0   1   1   4   8
   0   2   4   2   6
   0   2   4   2   5
   0   1   1   5   5
~~~
{:class="out"}

The index`(1:3:10, 1:2:10)` means "Rows 1 through 10 in steps of 
3 and columns 1 through 10 in steps of 2". Matlab will stop when it 
reaches or crosses the upper bounds of this index, i.e., it will
access columns 1, 3, 5, 7, and 9, but not 11.


<!-- keyword -->

The `:` by itself can be used to to slice an entire row or column. For
example, to get the fifth row of `patient_data`, we can do:

~~~
patient_data(5, :)
~~~
{:class="in"}


~~~
 Columns 1 through 30:

    0    1    1    3    3    1    3    5    2    4    4    7    6    5    3   10    8   10    6   17    9   14    9    7   13    9   12    6    7    7

 Columns 31 through 40:

    9    6    3    2    2    4    2    0    1    1
~~~
{:class="out"}


Finally, we can use the `end` keyword to refer to the end of a row
or column:

~~~
patient_data(50:end, 7)
~~~
{:class="in"}

~~~
ans =

   3
   1
   3
   4
   1
   1
   2
   2
   4
   4
   3
~~~
{:class="out"}

Matlab knows how to perform common mathematical operations on arrays.
If we want to find the average inflammation for all patients on all days,
we can just ask for the mean of the array:


~~~
mean(patient_data(:))
~~~
{:class="in"}

~~~
ans = 6.1487
~~~
{:class="out"}

The reason we couldn't just do `mean(patient_data)` is because, that 
would compute the mean of *each column* in our table, and return a list
of mean values. The expression `patient_data(:)` *flattens* the table intto a
one-dimensional array.

To get details about what a function, like `mean`,
does and how to use it, use Matlab's `help` command.


~~~
help mean
~~~
{:class="in"}

~~~
 -- Function File: mean (X)
 -- Function File: mean (X, DIM)
 -- Function File: mean (X, OPT)
 -- Function File: mean (X, DIM, OPT)
     Compute the mean of the elements of the vector X.

          mean (x) = SUM_i x(i) / N

     If X is a matrix, compute the mean for each column and return them
     in a row vector.

     The optional argument OPT selects the type of mean to compute.  The
     following options are recognized:

     "a"
          Compute the (ordinary) arithmetic mean.  [default]

     "g"
          Compute the geometric mean.

     "h"
          Compute the harmonic mean.

     If the optional argument DIM is given, operate along this
     dimension.

     Both DIM and OPT are optional.  If both are supplied, either may
     appear first.

     See also: median, mode.
~~~
{:class="out"}


~~~
disp(['Maximum inflammation: ', num2str(max(patient_data(:)))]);
disp(['Minimum inflammation: ', num2str(min(patient_data(:)))]);
disp(['Standard deviation: ', num2str(std(patient_data(:)))]);
~~~
{:class="in"}

~~~
Maximum inflammation: 20
Minimum inflammation: 0
Standard deviation: 4.6148
~~~
{:class="out"}

When analyzing data, though, we often want to look at partial statistics,
such as the maximum value per patient or the average value per day.
One way to do this is to assign the data we want to a new temporary 
array, then ask it to do the calculation:

~~~
patient_1 = patient_data(1, :)
disp(["Maximum inflation for patient 1: ", max(patient_1));
~~~
{:class="in"}

~~~
Maximum inflation for patient 1: 18
~~~
{:class="out"}

We don't actually need to store the row in a variable of its own.
Instead, we can combine the selection and the function call:

~~~
max(patient_data(1, :))
~~~
{:class="in"}

~~~
ans = 18
~~~
{:class="out"}


What if we need the maximum inflammation for *all* patients, or the
average for each day?
As the diagram below shows, we want to perform the operation across an
axis:

<div>
  <img src="img/matlab-operations-across-axes.svg" alt="Operations Across Axes" />
</div>

To support this, Matlab allows us to specify the *dimension* we 
want to work on. If we ask for the average across the dimension 1, 
we get:

~~~
mean(patient_data, 1)
~~~
{:class="in"}


~~~
ans =

 Columns 1 through 13:

    0.00000    0.45000    1.11667    1.75000    2.43333    3.15000    3.80000    3.88333    5.23333    5.51667    5.95000    5.90000    8.35000

 Columns 14 through 26:

    7.73333    8.36667    9.50000    9.58333   10.63333   11.56667   12.35000   13.25000   11.96667   11.03333   10.16667   10.00000    8.66667

 Columns 27 through 39:

    9.15000    7.25000    7.33333    6.58333    6.06667    5.95000    5.11667    3.60000    3.30000    3.56667    2.48333    1.50000    1.13333

 Column 40:

    0.56667

~~~
{:class="out"}


As a quick check, we can ask this array what its shape is:

~~~
size(mean(patient_data, 1))
~~~
{:class="in"}

~~~
ans = 
    1    40
~~~
{:class="out"}

The shape tells us we have a 1-by-40 vector, so this is the average
inflammation per day for all patients. If we average across axis 2, we
get:


~~~
mean(patient_data, 2)
~~~
{:class="in"}

~~~
ans =

   5.4500
   5.4250
   6.1000
   5.9000
   5.5500
   6.2250
   5.9750
   6.6500
   6.6250
   6.5250
   6.7750
   5.8000
   6.2250
   5.7500
   5.2250
   6.3000
   6.5500
   5.7000
   5.8500
   6.5500
   5.7750
   5.8250
   6.1750
   6.1000
   5.8000
   6.4250
   6.0500
   6.0250
   6.1750
   6.5500
   6.1750
   6.3500
   6.7250
   6.1250
   7.0750
   5.7250
   5.9250
   6.1500
   6.0750
   5.7500
   5.9750
   5.7250
   6.3000
   5.9000
   6.7500
   5.9250
   7.2250
   6.1500
   5.9500
   6.2750
   5.7000
   6.1000
   6.8250
   5.9750
   6.7250
   5.7000
   6.2500
   6.4000
   7.0500
   5.9000
~~~
{:class="out"}

which is the average inflammation per patient across
all days. 


### Plotting

The mathematician Richard Hamming once said,
"The purpose of computing is insight, not numbers," and the best
way to develop insight is often to visualize data. Visualization
deserves an entire lecture (or course) of its own, but we can 
explore a few features of Matlab here.

Let's display a heat map of our data:

~~~
imagesc(patient_data)
~~~
{:class="in"}

<div class="out">
<img src="img/01-intro_1.png" style="height:250px">
</div>

The `imagesc` function represents the matrix as a color image. Every
value in the matrix is *mapped* to a color. Blue regions in this heat map
are low values, while red shows high values.
As we can see,
inflammation rises and falls over a 40 day period.
Let's take a look at the average inflammation over time:

~~~
ave_inflammation = mean(patient_data, 1);
plot(ave_inflammation);
~~~
{:class="in"}

<div class="out">
<img src="img/01-intro_2.png" style="height:250px">
</div>

Here, we have put the average per day across all patients in the 
variable `ave_inflammation`, then used the `plot` function to display
a line graph of those values.
The result is roughly a linear rise and fall,
which is suspicious:
based on other studies, we expect a sharper rise and slower fall.
Let's have a look at two other statistics: the maximum and minimum
inflammation per day across all patients.

~~~
plot(max(data, [], 1));
title("Maximum inflammation per day");
plot(min(data, [], 1));
title("Minimum inflammation per day");
~~~
{:class="in"}

<div class="out">
<img src="img/01-intro_3.png" style="height:300px">
</div>

<div class="out">
<img src="img/01-intro_4.png" style="height:300px">
</div>

Like `mean()`, the functions
`max()` and `min()` can also operate across a specified dimension of
the matrix. However, the syntax is slightly different. To see why, 
run a `help` on each of these functions.

From the figures, we see that the maximum value rises and falls perfectly
smoothly, while the minimum seems to be a step function. Neither result
seems particularly likely, so either there 's a mistake in our 
calculations or something is wrong with our data.


#### Wrapping up

It's common to put multiple figures "side-by-side" in a single 
window for presentation and convenience. Here's how to use
the `subplot` function to do this:

~~~
subplot(1, 2, 1);
plot(max(patient_data, [], 1));
ylabel('max')

subplot(1, 2, 2);
plot(min(patient_data, [], 2));
ylabel('min')
~~~
{:class="in"}

<div class="out">
<img src="img/01-intro_5.png" style="width:750px; height:250px">
</div>


<div class="keypoints" markdown="1">
#### Key Points

* A value is a piece of data that a program works with. In Matlab,
values are always stored as arrays. Values can be numbers, vectors,
matrices, characters, strings, or a number of other things, but they are
*always* represented by an array.
* The "size" of an array is the number of elements it has in each
dimension.
* Use `variable = value` to assign a value to a variable. The "value"
is always an array in Matlab, and can contain different kinds 
of data, like integers, floats or characters.
* Use the `disp()` function to print things to screen. If you want
to print many things on the same line, print an array of things.
* Use `array(x, y)` to select a single element from an array.
* Array indices in Matlab start at 1, not 0 like in C or Python.
* Use `low:high` to specify a slice that includes the indices
from `low` to `high`.
* `mean()`, `max()` and `min()` can be used to calculate simple
statistics about an array, or along specified dimensions of that
array.
* The `imagesc()` function produces a heat map 
from a matrix (2-D array), and the `plot()` function can 
be used to produce a line graph from a vector (1-D array).

</div>

#### Next Steps
Our work so far has convinced us that something is wrong with our
first data file. We would like to check the other 11 the same way,
but typing in the same commands repeatedly is tedious and error-prone.
Since computers don't get bored (that we know of), we should create a
way to do a complete analysis with a single command, and then figure out
how to repeat that step once for each file. These operations are the 
subjects of the next two lessons.

