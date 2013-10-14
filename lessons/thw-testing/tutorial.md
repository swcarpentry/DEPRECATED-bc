---
layout: lesson
root: ../..
title: Testing Software
---
**Based on materials by Katy Huff, Rachel Slaybaugh, and Anthony Scopatz**

Introduction
------------

Now that you understand the basics of programming in Python, we'll move on to 
discuss two topics in "software engineering", which are how to test your code 
for accuracy and how to turn your code into stand alone scripts, or programs, 
that you can run from the command line.

Unit Testing Concepts
---------------------

As practicing scientists, we would never trust a lab measurement that we made 
with uncalibrated instruments. Similarly, as computational scientists, we 
shouldn't trust the results that our code gives us until we have tested it. 
Without calibration/testing, how do we know that our code is giving us the 
right answers?

In this lesson, we'll focus on unit tests, perhaps the most basic type of 
testing that we can run. Unit tests focus on a single "unit" of code, which in 
our case will be functions that we've written. We'll write tests to ensure that 
when our function is given a certain set of arguments as input, it generates 
output that we know to be correct. Once we have a complete test suite for a 
function, we can run the entire suite to make sure that all the tests pass (ie, 
that our function gives the correct output for all the combinations of input 
that we have decided to test).

For example, let's say that we have a function that reads in a data file, does 
some processing, and returns a result. We can test the function by giving it a 
small data file, for which we can calculate the correct result by hand, and 
making sure that the function gives the correct answer for this small file. 
This gives us more confidence that if we run the function on a different data 
set, perhaps a huge one for which we can't verify the results by hand, that 
we'll get an accurate result.

Even better, if we make changes to the internals of our function, we can run 
our tests again to make sure that we haven't accidentally broken anything (this 
is known as a "regression"). This makes us more free to continue to improve the 
performance of our code over time, and avoids the dreaded "it's working, don't 
touch it" phenomena.

In this lesson, we're going to use the simple and very popular `nose` package 
to write and run our tests.

A Unit Testing Example
----------------------

We'll practice unit testing using a function that we've already written to 
extract the mean number of animals seen per sighting from a csv file. First, 
let's place this function in an external module. To do this, copy the code 
below into a text file in this directory, and name it `mean_sightings.py`.

	import matplotlib.mlab as ml
	import numpy as np
	
	def get_sightings(filename, focusanimal):
	
		# Load table
		tab = ml.csv2rec(filename)
	
		# Find number of records and total count of animals seen
		isfocus = (tab['animal'] == focusanimal)
		totalrecs = np.sum(isfocus)
		meancount = np.mean(tab['count'][isfocus])
	
		# Return num of records and animals seen
		return totalrecs, meancount
	
This function uses boolean arrays to calculate the total number of records and 
mean number of animals per sighting for the focus animal.

To confirm that everything's working correctly, open up a new IPython notebook 
(in this same directory) and run the following in a cell:

	from mean_sightings import get_sightings
	print get_sightings('sightings_tab_sm.csv', 'Owl')

This should give you the correct answer for the Owl (check to make sure by 
looking at the raw data file and counting by hand).

Now that we have the function in a module, let's write some unit tests to make 
sure that the function is giving us the correct answers. Create a new text file 
called `test_mean_sightings.py`, which will hold our unit tests. At the top of 
this file, type (or copy) in the following code, which will import the function 
that we wish to test and set the filename that we want to use for the testing.

	from mean_sightings import get_sightings

	filename = 'sightings_tab_sm.csv'

Note that we are using a small, "toy" data set for testing so that we can 
calculate correct answers by hand.

Now, let's write our first test function, which will simply test to make sure 
that our function gives the correct answer when called using this small data 
set and the Owl as arguments. Test functions (written for the `nose` testing 
package) can contain any type of Python code, like regular functions, but have 
a few key features. First, they don't take any arguments. Second, they contain 
at least one `assert` statement - the test will pass if the condition following 
the `assert` statement is True, and the test will fail if it's False.

An example will make this more clear. Here's a test that checks whether the 
function returns the correct answers for the small data set and the Owl. Copy 
and paste this at the end of the `test_mean_sightings.py` file.

	def test_owl_is_correct():
	    owlrec, owlmean = get_sightings(filename, 'Owl')
		assert owlrec == 2, 'Number of records for owl is wrong'
	    assert owlmean == 17, 'Mean sightings for owl is wrong'

Note that we calculated the correct values of `owlrec` and `owlmean` by hand. 
Make sure that you get these right!

Now we're ready to run our suite of tests (so far, just this one test). Open a 
command line window, and `cd` to the directory containing your new Python 
files. Type `nosetests`, and examine the output. It should look something like 
this:

	.
	----------------------------------------------------------------------
	Ran 1 test in 0.160s
	
	OK

The dot on the first line shows that we had one test, and that it passed. There 
is one character printed for each test. A '.' means the test passed, a 'F' 
means the test failed, and an 'E' means there was an error in the test function 
itself.

Just for fun, try changing your test so that it fails (for example, assert that 
the number of Owl records should be 3). What output do you see now? Don't 
forget to change the test back so that it passes after you're done.

>### Exercise 1 - Test the Muskox results
>
>Add an additional test to your test file to make sure that your function also 
>gives the right answer when the animal is a Muskox. Run `nosetests` and make 
>sure both tests pass.

Great, now we have two tests that pass. However, both of these tests were 
fairly straightforward, in that they tested the expected behavior of the 
function under "normal" inputs. What about corner or boundary cases? For 
example, what should our function do if the animal is not found anywhere in the 
data set?

Let's say that we decide that our function should return 0 for the number of 
records and 0 for the mean animals per record if the animal is not found in the 
data set. Let's write a test to see if our function does this already:

	def test_animal_not_present():
	    animrec, animmean = get_sightings(filename, 'NotPresent')
		assert animrec == 0, 'Animal missing should return zero records'
	    assert animmean == 0, 'Animal missing should return zero mean'

If we run our test suite now, we see that this test fails. The output doesn't 
give us much of a hint as to what went wrong though - we know that animmean was 
not equal to zero, but what was it?

To find out, add the line `print animrec, animmean` right above the first 
assert statement, run the test suite again, and look at the output. Now we can 
see that the animmean was 'nan', which stands for "not a number". This is 
because when an animal is not found, our current function returns 0 for the 
number of records and 0 for the total count. To calculate the mean, it tries to 
divide 0/0, and gets 'nan'.

>### Exercise 2 - Fixing our function for a boundary case
>
>Modify the function `get_sightings` so that if the animal is not present, both 
>totalrecs and meancount are 0. HINT: Check if totalrecs is zero before 
>calculating meancount - if totalrecs is zero, meancount must also be zero.
>
>Run your test suite again to make sure all three tests now pass.

Here's another special case - all of the animal names in the data sets are 
capitalized, with the first letter in uppercase and the rest of the letters in 
lowercase. What if someone enters the name of the animal using the wrong case. 
For example, they might call the function with the argument 'oWl' for the 
animal name.

>### Exercise 3 - Fixing our function for bad input
>
>Write a test function that will pass only if your function returns the correct
>answer for owls if the input argument focusanimal is set to 'oWl'. Run this 
>test, and see that it currently fails.
>
>Then, modify the function so that this test passes. HINT: You can use the 
>method 'capitalize' on any string to correct its capitalization.
>
>Run your test suite again to make sure all four tests now pass.
>
>__Bonus__
>
>Determine what your function should return if a user gives the function a file 
>that does not exist. Write a test that checks that this value is indeed 
>returned for the case of a missing file, and modify your function to return it 
>as desired.

You can imagine adding more test functions as you think of more unusual cases 
that you want your function to correctly address. It is not unusual for the 
file containing test cases to be several times longer than the file containing 
the actual functions!

Now we're in a great position - we now have more confidence that our code is 
doing what we expect it to do.

Now let's say that we are planning to share our code with a colleague who is 
less experienced with programming, and we think that he/she might not 
understand the neat boolean indexing tricks that we've been using. For clarity, 
we decide that we'll replace the guts of our `get_sightings` function with code 
that calculates the same thing but uses a for loop instead. We've already 
written this code in the previous lesson, so we can simply erase our existing 
`get_sightings` function and replace it with this code instead:


	def get_sightings(filename, focusanimal):
	
	    # Load table
	    tab = ml.csv2rec(filename)
		
		# Standardize capitalization of focusanimal
		focusanimal = focusanimal.capitalize()

	    # Loop through all records, countings recs and animals
	    totalrecs = 0
	    totalcount = 0
	    for rec in tab:
			if rec['animal'] == focusanimal:
	            totalrecs += 1
	            totalcount += rec['count']
	
	    meancount = totalcount/totalrecs
	
	    # Return num of records and animals seen
	    return totalrecs, meancount

Thinking ahead, we made sure to add a line to fix the capitalization problem 
right away so that our fourth unit test should pass. Since this code worked 
before, we're confident that it will work now. Just to be sure, though we run 
our test suite again.

>### Exercise 4 - Examining and fixing regressions
>
>You are shocked to discover that two of the four tests now fail! How can this 
>be? We were sure that the new for loop code was correct, and we looked at its 
>output before to convince ourselves that it was correct...
>
>Try to uncover the causes of this regression. One failure should have a fairly 
>obvious cause (it relates to the issue of an animal not being present, which 
>we check with the third test). The second failure has a more subtle cause - 
>try to figure out the problem, and correct the function to give the right 
>answer.

### Test Driven Development - the joy of Red/Green/Refactor

Instead of fixing the above code, we're going to delete get_sightings, and do a very simple run through TDD.

The big idea here is that you think about your problem and write your unit tests *before* 
you write a single line of code. 
- This forces you to think about what your problem in terms of different modes of 
  success/failure and various edge cases, rather than just the basic functionality. 
- It means that you implement the right amount of functionality without overbuilding.
- It also gives you a ready-made specification for your design

We have already written our first 4 test cases. 
- Run ``nosetests``. You will see everything fail (Red)

Now we're going to write a bare minimum ``get_sightings`` that passes the first test case. The code will be 
really stupid

    def get_sightings(filename, focusanimal):
		return (2, 17)
	
This is clearly wrong BUT it passes a couple of test cases. It has also forced you to think about the structure of your function. 

Now that you have a couple of Greens you would refactor the code to be a little smarter. 

Continue to repeat this process of turning Red to Green; then refactoring and cleaning up.

Hopefully, this actually helps you write better code that has fewer bugs, and gives you deeper insight into the structure of your 
program.

Example:

    def get_sightings(filename, focusanimal):
    
    	# Load table
    	tab = ml.csv2rec(filename)
    
    	# Standardize capitalization of focusanimal
    	focusanimal = focusanimal.capitalize()
    
    	# Loop through all records, countings recs and animals
        totalrecs = 0.
        totalcount = 0.
    	for rec in tab:
            if rec['animal'] == focusanimal:
            	totalrecs += 1
            	totalcount += rec['count']
    
    	if totalrecs==0:
            meancount = 0
    	else:
        	meancount = totalcount/totalrecs
    
    	# Return num of records and animals seen
    	return totalrecs, meancount

__BONUS__ If there is time, write some tests that will pass for a different csv file.

Making a Standalone Script
--------------------------

Now that our module has been tested, lets turn this program into a standalone 
script that we can run from the command line. This takes very little additional 
work, now that we have our function in a module.

At the bottom of the `mean_sightings.py`, add the following lines:

	filename = 'sightings_tab_sm.csv'
	focusanimal = 'Owl'
	print get_sightings(filename, focusanimal)

Now, head over to the command line and make sure that you're in the directory 
containing the `mean_sightings.py` file. Type the statement below then hit 
return.

	python mean_sightings.py

You should see the output `(2, 17)` printed to the screen, which is the correct 
number of records and the mean number of animals per record for the Owl in the 
`sightings_tab_sm.csv` file.

This is interesting, but it would be much more useful if we could give our 
command line program arguments, in the same way that we would type `cat 
myfile.txt`. For example, we may want to type `python mean_sightings.py 
sightings_tab_sm.csv Owl` instead of having to make a change in the file itself 
each time we want to use a different file and focal animal.

This is actually pretty easy to do using a Python module called `sys`. At the 
top of the `mean_sightings.py` file, add the line

	import sys

then at the bottom of the file, change your code to read

	filename = sys.argv[1]
	focusanimal = sys.argv[2]
	print get_sightings(filename, focusanimal)

The variable `sys.argv` is a list of all of the arguments given on the command 
line when this file is called (you can see this by putting `print sys.argv` a 
the bottom of the script as well. The first argument, `sys.argv[0]`, is always 
the name of the file that was run - in this case, it's `mean_sightings.py`. The 
second and third arguments are stored in `sys.argv[1]` and `sys.argv[2]`, and 
we've chosen to use these as the filename and focusanimal.

Now you can simply type

	python mean_sightings.py sightings_tab_sm.csv Owl

and you'll get what you were expecting. Try this out with different animals and 
with the large table. Make sure it works for our special cases that we 
addressed before, like the capitalization of the animal name being incorrect.

Two more small changes will make our command line script extra professional.

First, we have now changed our file `mean_sightings.py` so that it runs from 
the command line, but what if we want to also be able to import functions from 
it as a module from other Python programs (such as in notebooks when we run 
`import mean_sightings`)? The best way to do this is to wrap all of the lines 
at the bottom of our file (the ones that produce the command line output, not 
the functions themselves) into a special if statement like so:

	if __name__ == '__main__':
		filename = sys.argv[1]
		focusanimal = sys.argv[2]
		print get_sightings(filename, focusanimal)

When a Python script is run from the command line, a special hidden variable 
called `__name__` is set to equal the string `__main__`. This special if 
statement thus encloses code that we only want to run when the file is run from 
the command line, not when it's imported by another file. You'll see this 
special statement in many Python scripts.

Second, we can set up our file so that it can be executed directly like any 
other shell script (so that we can run `mean_sightings` from the command line 
instead of `python mean_sightings`). To do this, we have to first tell our 
shell that when this file is executed directly, it should be run using the 
python interpreter. To do this, make the very first line of the file

	#!/usr/bin/env python

Then, we need to give the file `mean_animals.py` permission to execute on its 
own. From the command line, in the directory containing the file 
`mean_animals.py`, run the line

	chmod 755 mean_sightings.py

Now we can run our file as a standalone script simply by executing the 
statement

	./mean_sightings.py sightings_tab_sm.csv Owl

That annoying little `./` at the front is because the shell, by default, 
doesn't look inside your current directory for executable programs - it only 
looks for executables within directories specified by the PATH shell variable. 
We'll leave it as an exercise to you to look up how to add a directory to your 
PATH. If you have certain scripts that you run very often, a common trick is to 
create a single directory, such as `~/bin/`, add this to your PATH permanently 
by modifying your `.bashrc` or `.bash_profile`, and put all of your commonly 
executed scripts in that directory. Then you will be able to call them from the 
command line as `scriptname`, just like all of the built in shell commands.

As a side note, you are free to remove the extension from your script file 
names if you'd like. For example, you are free to rename `mean_sightings.py` to 
`mean_sightings` - everything will still work as expected.
