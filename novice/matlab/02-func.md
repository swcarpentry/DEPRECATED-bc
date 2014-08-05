---
layout: lesson
root: ../..
---


## Creating Functions

If we only had one data set to analyze,
it would probably be faster to load the file into a spreadsheet
and use that to plot some simple statistics.
But we have twelve files to check,
and may have more in future.
In this lesson,
we'll learn how to write a function
so that we can repeat several operations with a single command.

<div class="objectives">
### Objectives
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

A Matlab function *must* be saved in a text file with a `.m` extension. The name of that file must be the same as the function defined
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
writing a function named
`fahr_to_kelvin`, that accepts a single parameter, `ftemp`, and outputs a 
single value, `ktemp`. 

Anything following the function definition line is called the *body* of the
function.

We can call our function from the command line like any other function:
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


