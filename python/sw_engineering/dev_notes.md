Dev Notes
=========

Things we want to demonstrate:

- write lots of small tools (functions)
- document your code
- test your code, including TDD
- debugging
  -   build in a crash somewhere and explain tracebacks
  -   do an import pdb; pdb.set_trace()

Create a few utilities to work with the `*animal.txt` files from
introductory Python.

Put all the functionality in one file and make small scripts that import
from that file, parse command line arguments, and do stuff.

- average number of an animal seen per sighting

Students may want an IPython notebook open as a scratch pad.

Lesson Plan
-----------

We're going to make a command line script that will print the average
group size in a given file for a given animal.

0.  Discuss libraries and re-using code.
1.  Demonstrate importing a module (maybe `math` or `glob`), then
    demonstrate putting a function in a file and importing it.

    Exercise: Make a new text file called `animals.py`. Copy the file
    reading function from yesterday's IPython notebook into the file and
    modify it so that it returns the columns of the file as lists (instead
    of printing certain lines). (They may want to develop the function in
    the IPython notebook and then copy it over.)

2.  Go over the exercise and talk about documentation as you do.
3.  How do we know the function works correctly? Try importing and
    running it on a small file.
4.  But what if we want to make changes and make sure it still works
    afterward, or we want to make sure it isn't broken when we add new
    stuff later?
    - Demonstrate a potential test solution comparing output in an `if` and
      printing a message if something doesn't match.
    - Explain `assert` and demonstrate
    - Explain unit tests and show how to run with nosetests.

5.  Explain test driven development.
    Exercise: We're going to make a function to calculate the mean of all
    the values in a list, but we're going to write the tests for it first.
    Make a new text file called `test_animals.py`. Make a function called
    `test_mean` that runs your theoretical mean function through several
    tests.

6.  When going over this with the students, put in a test with an empty
    list. Also put in tests with lists that contain all ints.

    Exercise: Write the mean function in `animals.py` and verify that it
    passes your tests.

7.  When going over this with the students do not put in a test for an
    empty list. The error when running the tests will give a chance to
    teach tracebacks. Also do not put in any coercion to float so some
    means will be wrong due to integer truncation. More debugging,
    `--pdb-failure`.

8.  The last piece we'll need is a function that takes the output of the
    file reader function and returns only the data relevant to a given
    animal. Write this function as a live exercise with the students
    participating, though they can do it on their own if they want.

    Exercise: Write tests for a function that will take a file name and
    animal name as arguments, and return the average number of animals per
    sighting.

    Exercise: Write a function that takes a file name and animal name and
    returns the average number of animals per sighting. Make sure it passes
    your tests.

9.  After going over exercises, conclude by showing how to make a
    command line script that takes a file name and animal name as
    arguments and prints the average number of animals per sighting.

10. If time permits, could demonstrate `pdb.set_trace()` somewhere in
    the script execution.
