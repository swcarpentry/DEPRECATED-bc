---
layout: lesson
root: ../..
title: Exceptions
---
Assertions help us catch errors in our code,
but things can go wrong for other reasons,
like missing or badly-formatted files.
Most modern programming languages allow programmers to use
[exceptions](../../gloss.html#exception) to separate
what the program should do if everything goes right
from what it should do if something goes wrong.
Doing this makes both cases easier to read and understand.

For example,
here's a small piece of code that tries to read parameters and a grid from two
separate files,
and reports an error if either goes wrong:

~~~
try:
    params = read_params(param_file)
    grid = read_grid(grid_file)
except:
    log.error('Failed to read input file(s)')
    sys.exit(ERROR)
~~~
{:class="in"}

We join the normal case and the error-handling code using the keywords `try` and
`except`.
These work together like `if` and `else`:
the statements under the `try` are what should happen if everything works,
while the statements under `except` are what the program should do if something
goes wrong.

We have actually seen exceptions before without knowing it,
since by default,
when an exception occurs,
Python prints it out and halts our program.
For example,
trying to open a nonexistent file triggers a type of exception called an
`IOError`,
while trying to access a list element that doesn't exist
causes an `IndexError`:

~~~
open('nonexistent-file.txt', 'r')
~~~
{:class="in"}
~~~
---------------------------------------------------------------------------
IOError                                   Traceback (most recent call last)

<ipython-input-13-58cbde3dd63c> in <module>()
----> 1 open('nonexistent-file.txt', 'r')

IOError: [Errno 2] No such file or directory: 'nonexistent-file.txt'
~~~
{:class="err"}
~~~
values = [0, 1, 2]
print values[999]
~~~
{:class="in"}
~~~
---------------------------------------------------------------------------
IndexError                                Traceback (most recent call last)

<ipython-input-14-7fed13afc650> in <module>()
1 values = [0, 1, 2]
----> 2 print values[999]

IndexError: list index out of range
~~~
{:class="err"}

We can use `try` and `except` to deal with these errors ourselves
if we don't want the program simply to fall over:

~~~
try:
    reader = open('nonexistent-file.txt', 'r')
except IOError:
    print 'Whoops!'
~~~
{:class="in"}
~~~
Whoops!
~~~
{:class="err"}

When Python executes this code,
it runs the statement inside the `try`.
If that works, it skips over the `except` block without running it.
If an exception occurs inside the `try` block,
though,
Python compares the type of the exception to the type specified by the `except`.
If they match, it executes the code in the `except` block.

`IOError` is the particular kind of exception Python uses
when there is a problem related to input and output,
such as files not existing
or the program not having the permissions it needs to read them.
We can put as many lines of code in a `try` block as we want,
just as we can put many statements under an `if`.
We can also handle several different kinds of errors afterward.
For example,
here's some code to calculate the entropy at each point in a grid:

~~~
try:
    params = read_params(param_file)
    grid = read_grid(grid_file)
    entropy = lee_entropy(params, grid)
    write_entropy(entropy_file, entropy)
except IOError:
    report_error_and_exit('IO error')
except ArithmeticError:
    report_error_and_exit('Arithmetic error')
~~~
{:class="in"}

Python tries to run the four functions inside the `try` as normal.
If an error occurs in any of them,
Python immediately jumps down
and tries to find an `except` of the corresponding type:
if the exception is an `IOError`,
Python jumps into the first error handler,
while if it's an `ArithmeticError`,
Python jumps into the second handler instead.
It will only execute one of these,
just as it will only execute one branch
of a series of `if`/`elif`/`else` statements.

This layout has made the code easier to read,
but we've lost something important:
the message printed out by the `IOError` branch doesn't tell us
which file caused the problem.
We can do better if we capture and hang on to the object that Python creates
to record information about the error:

~~~
try:
    params = read_params(param_file)
    grid = read_grid(grid_file)
    entropy = lee_entropy(params, grid)
    write_entropy(entropy_file, entropy)
except IOError as err:
    report_error_and_exit('Cannot read/write' + err.filename)
except ArithmeticError as err:
    report_error_and_exit(err.message)
~~~
{:class="in"}

If something goes wrong in the `try`,
Python creates an exception object,
fills it with information,
and assigns it to the variable `err`.
(There's nothing special about this variable name&mdash;we can use anything we
want.)
Exactly what information is recorded depends on what kind of error occurred;
Python's documentation describes the properties of each type of error in detail,
but we can always just print the exception object.
In the case of an I/O error,
we print out the name of the file that caused the problem.
And in the case of an arithmetic error,
printing out the message embedded in the exception object is what Python would
have done anyway.

So much for how exceptions work:
how should they be used?
Some programmers use `try` and `except` to give their programs default
behaviors.
For example,
if this code can't read the grid file that the user has asked for,
it creates a default grid instead:

~~~
try:
    grid = read_grid(grid_file)
except IOError:
    grid = default_grid()
~~~
{:class="in"}

Other programmers would explicitly test for the grid file,
and use `if` and `else` for control flow:

~~~
if file_exists(grid_file):
    grid = read_grid(grid_file)
else:
    grid = default_grid()
~~~
{:class="in"}

It's mostly a matter of taste,
but we prefer the second style.
As a rule,
exceptions should only be used to handle exceptional cases.
If the program knows how to fall back to a default grid,
that's not an unexpected event.
Using `if` and `else`
instead of `try` and `except`
sends different signals to anyone reading our code,
even if they do the same thing.

Novices often ask another question about exception handling style,
but before we address it,
there's something in our example that you might not have noticed.
Exceptions can actually be thrown a long way:
they don't have to be handled immediately.
Take another look at this code:

~~~
try:
    params = read_params(param_file)
    grid = read_grid(grid_file)
    entropy = lee_entropy(params, grid)
    write_entropy(entropy_file, entropy)
except IOError as err:
    report_error_and_exit('Cannot read/write' + err.filename)
except ArithmeticError as err:
    report_error_and_exit(err.message)
~~~
{:class="in"}

The four lines in the `try` block are all function calls.
They might catch and handle exceptions themselves,
but if an exception occurs in one of them that *isn't* handled internally,
Python looks in the calling code for a matching `except`.
If it doesn't find one there,
it looks in that function's caller,
and so on.
If we get all the way back to the main program without finding an exception
handler,
Python's default behavior is to print an error message like the ones we've been
seeing all along.

This rule is the origin of the rule
[throw low, catch high](../rules.html#throw-low-catch-high).
There are many places in our program where an error might occur.
There are only a few, though, where errors can sensibly be handled.
For example,
a linear algebra library doesn't know whether it's being called directly from
the Python interpreter,
or whether it's being used as a component in a larger program.
In the latter case,
the library doesn't know if the program that's calling it is being run from the
command line or from a GUI.
The library therefore shouldn't try to handle or report errors itself,
because it has no way of knowing what the right way to do this is.
It should instead just [raise](../../gloss.html#raise) an exception,
and let its caller figure out how best to handle it.

Finally,
we can raise exceptions ourselves if we want to.
In fact,
we *should* do this,
since it's the standard way in Python to signal that something has gone wrong.
Here,
for example,
is a function that reads a grid and checks its consistency:

~~~
def read_grid(grid_file):
    data = read_raw_data(grid_file)
    if not grid_consistent(data):
        raise Exception('Inconsistent grid: ' + grid_file)
    result = normalize_grid(data)
    return result
~~~
{:class="in"}

The `raise` statement creates a new exception with a meaningful error message.
Since `read_grid` itself doesn't contain a `try`/`except` block,
this exception will always be thrown up and out of the function,
to be caught and handled by whoever is calling `read_grid`.
We can define new types of exceptions if we want to.
And we should,
so that errors in our code can be distinguished from errors in other people's
code.
However,
this involves classes and objects,
which is outside the scope of these lessons.
