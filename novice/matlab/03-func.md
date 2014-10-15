---
layout: lesson
root: ../..
title: Creating Functions
---

If we only had one data set to analyze,
it would probably be faster to load the file into a spreadsheet
and use that to plot some simple statistics.
But we have twelve files to check,
and may have more in future.
In this lesson,
we'll learn how to write a function
so that we can repeat several operations with a single command.

<div class="objectives" markdown="1">
#### Objectives
* Explain a Matlab function file
* Define a function that takes parameters
* Test and debug a function.
* Explain what a call stack is, and trace changes to the call stack as
    functions ar ecalled.
* Set default values for function parameters.
* Explain why we should divide programs into small, single-purpose functions.
</div>

### Defining a Function

Let's start by defining a function `fahr_to_kelvin` that converts temperatures from Fahrenheit to Kelvin:

~~~
% file fahr_to_kelvin.m

function ktemp = fahr_to_kelvin(ftemp)
    ktemp = ((ftemp - 32) * (5/9)) + 273.15;
end
~~~
{:class="in"}

A Matlab function *must* be saved in a text file with a `.m` extension.
The name of that file must be the same as the function defined
inside it. So, you will need to save the above code in a file called
`fahr_to_kelvin.m`.

<!-- FIXME: 00-basics lesson should talk about .m files, and
    the variable `ans` -->

<!-- FIXME: Nothing about Comments? -->

<!-- FIXME: this lesson should talk about multiple functions in a
    single file -->

The first line of our function:

~~~
function ktemp = fahr_to_kelvin(ftemp)
~~~

is called the *function definition*, and it declares that we're
writing a function named `fahr_to_kelvin`, that accepts a single parameter,
`ftemp`, and outputs a  single value, `ktemp`.  Anything following the
function definition line is called the *body* of the
function. The keyword `end` marks the end of the function body, and the
function won't know about any code after `end`.

Just as we saw with scripts, functions must be _visible_ to MATLAB, i.e.,
a file containing a function has to be placed in a directory that
MATLAB knows  about. The most convenient of those directories is the
current working directory.

We can call our function from the command line
like any other MATLAB function:

~~~
fahr_to_kelvin(32)
~~~
{:class="in"}

~~~
ans = 273.15
~~~
{:class="out"}

When we pass a value, like `32`, to the function, the value is assigned
to the variable `ftemp` so that it can be used inside the function. If we
want to return a value from the function, we must assign that value to a
variable named `ktemp`---in the first line of our function, we promised
that the output of our function would be named `ktemp`.

Outside of the function, the names `ftemp` and `ktemp` don't matter,
they are only used by the function body to refer to the input and
output values.


<!-- FIXME: make up a debugging scenario -->

### Composing Functions

Now that we've seen how to turn Fahrenheit to Kelvin, it's easy to turn
Kelvin to Celsius.

~~~
% file kelvin_to_celsius.m

function ctemp = kelvin_to_celsius(ktemp)
    ctemp = ktemp - 273.15;
end
~~~
{:class="in"}

Again, we can call this function like any other:

~~~
kelvin_to_celsius(0.0)
~~~
{:class="in"}

~~~
ans = -273.15
~~~
{:class="out"}

What about converting Fahrenheit to Celsius?
We could write out the formula, but we don't need to.
Instead, we can [compose](../../gloss.html#function-composition) the two
functions we have already created:

~~~
% file fahr_to_celsius.m

function ctemp = fahr_to_celsius(ftemp)
    ktemp = fahr_to_kelvin(ftemp);
    ctemp = kelvin_to_celsius(ktemp);
end
~~~
{:class="in"}

Calling this function,

~~~
kelvin_to_celsius(0.0)
~~~

we get, as expected:

~~~
ans = -273.15
~~~

This is our first taste of how larger programs are built:
we define basic operations,
then combine them in ever-large chunks to get the effect we want.
Real-life functions will usually be larger than the ones shown
here---typically half a dozen to a few dozen lines---but
they shouldn't ever be much longer than that,
or the next person who reads it won't be able to understand what's going on.

#### Challenges 

1. In Matlab, we concatenate strings by putting them into an array or using the
   `strcat` function:

   ~~~
   disp(['abra', 'cad', 'abra'])
   ~~~
   {:class="in"}
   ~~~
   abracadabra 
   ~~~
   {:class="out"}

   ~~~
   disp(strcat('a', 'b'))
   ~~~
   {:class="in"}
   ~~~
   ab
   ~~~
   {:class="out"}

   Write a function called `fence` that takes two parameters, `original` and
   `wrapper` and appends `wrapper` before and after `original`:

   ~~~
   disp(fence('name', '*'))
   ~~~
   {:class="in"}
   ~~~
   *name*
   ~~~
   {:class="out"}

1. If the variable `s` refers to a string, then `s(1)` is the string's first
   character and `s(end)` is its last. Write a function called `outer` that returns
   a string made up of just the first and last characters of its input: 

   ~~~
   disp(outer('helium'))
   ~~~
   {:class="in"}
   ~~~
   hm
   ~~~
   {:class="out"}

### The Call Stack


Let's take a closer look at what happens when we call
`fahr_to_celcius(32.0)`.
To make things clearer, we'll start by putting the initial value 32.0
in a variable and store the final result in one as well:

~~~
original = 32.0;
final = fahr_to_celcius(original);
~~~
{:class="in"}

The diagram below shows what memory looks like after the
first line has been executed:

<!-- FIXME: Writeup about Call Stacks -->



### Testing and Documentation

Once we start putting things in functions so that we can
re-use them, we need to start testing that those functions are
working correctly.
To see how to do this, let's write a function to center a
dataset around a particular value:

~~~
function out = center(data, desired)
    out = (data - mean(data)) + desired
end
~~~
{:class="in"}

We could test this on our actual data, but since we
don't know what the values ought to be,
it will be hard to tell if the result was correct,
Instead, let's use to create a matrix of 0's, and then center that
around 3:

~~~
z = zeros(2,2);
center(z, 3)
~~~
{:class="in"}

~~~
ans =

   3   3
   3   3
~~~
{:class="out"}

That looks right, let's try center on our real data:

~~~
data = csvread('inflammation-01.csv');
centered = center(data(:), 0)
~~~
{:class="in"}

It's hard to tell from the default output whether the
result is correct--this is often the case when working with
fairly large arrays--but, there are a few simple tests that
will reassure us.

Let's calculate some simple statistics:

~~~
disp([min(data(:)), mean(data(:)), max(data(:))])
~~~
{:class="in"}

~~~
0.00000    6.14875   20.00000
~~~
{:class="out"}

And let's do the same after applying our `center` function
to the data:

~~~
disp([min(centered(:)), mean(centered(:)), max(centered(:))])
~~~
{:class="in"}

~~~
-6.1487e+00  -2.2962e-14   1.3851e+01
~~~
{:class="out"}

<!-- FIXME: challenge -->

That seems almost right, the original mean
was about 6.1, so the lower bound from zero is now about -6.1.
The mean of the centered data isn't quite zero--we'll explore
why not in the challenges--but it's pretty close. We can even
go further and check that the standard
deviation hasn't changed:


~~~
std(data(:)) - std(centered)
~~~
{:class="in"}

~~~
5.3291e-15
~~~
{:class="out"}


The difference is very small. It's still possible that our function
is wrong, but it seems unlikely enough that we should probable
get back to doing our analysis. We have one more task first, though:
we should write some [documentation](../../gloss.html#documentation)
for our function to remind ourselves later what it's for and
how to use it.

~~~
function out = center(data, desired)
    %   Center data around a desired value.
    %
    %       center(DATA, DESIRED)
    %
    %   Returns a new array containing the values in
    %   DATA centered around the value.

    out = (data  - mean(data)) + desired;
end
~~~
{:class="in"}

Comment lines immediately below the function definition line
are called "help text". Typing `help function_name` brings
up the help text for that function:

~~~
help center
~~~
{:class="in"}

~~~
Center data around a desired value.

    center(DATA, DESIRED)

Returns a new array containing the values in
DATA centered around the value.
~~~
{:class="out"}

#### Challenges

1. Write a function called `rescale` that takes an array as input and returns an
   array of the same shape with its values scaled to lie in the range 0.0 to 1.0.
   (If L and H are the lowest and highest values in the input array, respectively,
   then the function should map a value v to (v - L)/(H - L).) Be sure to give the
   function a comment block explaining its use.

1. Run `help linspace` to see how to use this function to generate
   regularly-spaced values. Use arrays like this to test your `rescale` function.

1. Write a function `run_analysis` that accepts a filename
   as parameter, and displays the three graphs produced in the
   previous lesson, i.e., `run_analysis('inflammation-01.csv')`
   should produce the corresponding graphs for the first
   data set. Be sure to give your function help text.

<div class="keypoints" markdown="1">
#### Key Points

* Define a function in a file called `function_name.m`.
* The first line of a function is the definition line and
is most generally written as
`function [out1, out2 ..] = function_name(in1, in2 ...)`
* The body of a function must end with the keyword `end`.
* Use `help function_name` to get information on how to use a function.
* Unbroken comment lines below the function definition line are displayed
as help text.

</div>

### Next Steps

<!--FIXME: need next steps here -->

