---
layout: lesson
root: ../..
---

## Analyzing Patient Data


We are studying inflammation in patients who have been given a new treatment for arthritis,
and need to analyze the first dozen data sets.
The data sets are stored in <a href="../../gloss.html#csv">comma-separated values</a> (CSV) format:
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


<!-- FIXME: add gloss entries for terms like "array" -->


<div class="objectives" markdown="1">
* Learn about Matlab arrays.
* Read tabular data from a file into a program.
* Assign values to variables.
* Select individual values and subsections from data.
* Perform operations on arrays of data.
* Display simple graphs.
</div>

### Loading Data

<!-- FIXME: high level -->

I/O, or reading and writing data, is essential to scientific computing, 
and admittedly, something that we'd rather not spend a lot of time 
thinking about. Matlab spares us the work by providing a number 
of high-level functions that accomplish these 
tasks efficiently, while hiding the grisly details.

To fetch the data from our CSV file into Matlab, type following 
command into the Matlab shell, and press Enter:

~~~
csvread('inflammation-01.csv')
~~~
{:class="in"}

To suppress the
output, simply put a semicolon at the end of your commnd:

~~~
csvread('inflammation-01.csv');
~~~
{:class="in"}

<!-- FIXME: gloss function call, parameters, -->

The expression `csvread(...)` is a 
[function call](../../gloss.html#function-call). 
Functions generally need [parameters](../../gloss.html#parameter)
to run.
In the case of the `csvread` function, we need to provide a single
parameter: the name of the file we want to read data from. This 
parameter needs to be a character string or 
[string](../../gloss.html#string), so we put it in quotes.

<!-- FIXME: below not entirely true, maybe footnote required? -->

Our call to `csvread` read our file, and printed the data inside 
to the screen. But we still have no way to modify those values
or compute with them. To do that, we need to assign the array to a
[variable](../../gloss.html#variable).

~~~
a = csvread('inflammation-01.csv');
~~~
{:class="in"}

A variable is just a name for a piece of data or *value*. 
Variable names must begin with a letter, and can contain 
numbers or underscores. Examples of valid variable names are 
`x`, `current_temperature` and `subject_id`. 

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


div>
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


v<div>
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
`csvread` and save its result"

~~~
data = csvread('inflammation-01.csv');
~~~
{:class="in"}

Matlab provides a function
to list all variables that have been assigned data:

~~~
who
~~~
{:class="in"}

~~~
Variables in the current scope:


data
~~~
{:class="out"}


<div class="challenges">
#### Challenges

* Draw diagrams showing what varuables refer to what values 
after each statement in the following program:

~~~
mass = 47.5;
age = 122;
mass = mass*2;
age = age - 20;
~~~

</div>


### Manipulating Data






