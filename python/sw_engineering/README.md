Building a Library of Code you Trust
====================================

Suppose we're going to be dealing a lot with these animal count files,
and doing many different kinds of analysis with them. In the
introduction to Python lesson we wrote a function that reads these files
but it's stuck off in an IPython notebook. We could copy and paste it
into a new notebook every time we want to use it but that gets tedious
and makes it difficult to add features to the function. The ideal
solution would be to keep the function in one spot and use it over and
over again from many different places. Python modules to the rescue!

We're going to move beyond the IPython notebook. Most Python code is
stored in `.py` files and then used in other `.py` files where it
has been pulled in using an `import` statement. Today we'll show you
how to do that.

Exercises
=========

Exercise 1
----------

Make a new text file called `animals.py`. Copy the file reading
function from yesterday's IPython notebook into the file and modify it
so that it returns the columns of the file as lists (instead of printing
certain lines).

Exercise 2
----------

We're going to make a function to calculate the mean of all the values
in a list, but we're going to write the tests for it first. Make a new
text file called `test\_animals.py`. Make a function called
`test\_mean` that runs your theoretical mean function through several
tests.

Exercise 3
----------

Write the mean function in `animals.py` and verify that it passes your
tests.

Exercise 4
----------

Write tests for a function that will take a file name and animal name as
arguments, and return the average number of animals per sighting.

Exercise 5
----------

Write a function that takes a file name and animal name and returns the
average number of animals per sighting. Make sure it passes your tests.
