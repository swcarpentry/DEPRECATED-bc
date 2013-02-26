# Debugging

**Materials contributed by Anthony Scopatz and Patrick Fuller**

## Exercise: What is debugging?

Before I show you the practice (ahem, _art_) of debugging, separate out into groups of 2-3 
people. Follow these steps:

1.  Come up with a definition of debugging.
2.  Write it down on a strip of paper.
3.  Give me the strip of paper.
4.  ???
5.  Profit.

(Bonus Challenge: Make a new friend!)

Time limit: 5 min.

## Why does debugging matter?

Unless you're perfect, you are bound to make errors. Especially early 
on, expect to spend much more time debugging than actually coding. The process 
fits the Pareto principle - you're going to spend \~20% of your time writing 
~80% of your code, and the other ~80% of your time will be spent screaming 
obscenities at your computer (I think that's what the Pareto principle says, 
anyway). Remember to keep calm, and **LEARN** from your mistakes.

## Debugging 101: exceptions, errors, and tracebacks

When your code errors, Python will stop and return an _exception_ that attempts 
to tell you what's up. There are \~165 exceptions in the Python standard 
library, and you'll be seeing many of them very soon. Exceptions to befriend 
include...

```python
SyntaxError # You're probably missing a parenthesis or colon
NameError   # There's probably a variable name typo somewhere
TypeError   # You're doing something with incompatible variable types
ValueError  # You're calling a function with the wrong parameter
IOError     # You're trying to use a file that doesn't exist
IndexError  # You're trying to reference a list element that doesn't exist
KeyError    # Similar to an IndexError, but for dictionaries
Exception   # This means "an error of any type" - hopefully you don't see it often
```

When code returns an exception, we say that the exception was _thrown_ or
_raised_. These exceptions may be _handled_ or _caught_ by the code. Speaking
of, you can handle exceptions in Python like so:

```python
try:
    a = 1.0 / 0.0
except ZeroDivisionError:
    print "Going from zero to hero."
    a = 1.0
```

That being said, there are some things you should keep in mind.
 * Exception handling in your own code should be seen as a last resort. _Never_
   use exception handling where another approach would work just as well.
 * If you have to handle exceptions, be specific in their type. Writing a blanket
   `except Exception` line provides a place for unintended bugs to hide.

When an exception is printed, it often comes with something called a
_traceback_. This is Python's attempt to tell you where the code errored. It
will look like gibberish for a while, but that impression will go away with time.

So, when your code errors, Python tells you 1. why it errored and 2. where it
errored. "Isn't that enough to debug?", you might ask. Well, yeah. It is. But,
if you debug only by running your code, you're going to be spending a lot more
time in the screaming-obscenities-at-your-computer portion of coding. Every tool
discussed below doesn't add much in terms of functionality (they're still just
pointing out errors), but they all help in decreasing debugging time.

## Linting: catching the stupid errors

As I said before, you can debug by simply attempting to run your code. This,
however, is very annoying. First off, the code will always stop at the first 
exception. This means that, if you have ten errors, you'll have to run the code 
ten times to find all of them. Now, imagine that this is long-running code.
Imagine waiting five minutes for your code to run, only to discover it breaks
because of a typo. Doesn't that sound terrible?

Enter linting. "Linting" is the process of discovering errors in a code (typically 
typos and syntax errors... ie. the dumb stuff) before the code is ever run or 
compiled. In Python, this can be accomplished through using the [pyflakes](http://pypi.python.org/pypi/pyflakes/) 
library. It works by statically analyzing your code without running it. This 
means that it can find multiple errors at once, rather than stopping at the 
first exception.

You can run pyflakes on your code by typing

```
pyflakes my_code.py
```

We can take this a step further with integrated development environments, or
_IDEs_. IDEs are (basically) glorified text editors that dynamically lint,
showing you typos as you write your code. You will find that coders generally
have strong opinions on the use of IDEs, either positive or negative. Regardless,
if you want to play around with one, I recommend [Eclipse](http://www.eclipse.org/) 
with the [PyDev](http://pydev.org/) plugin.

## Coding standards: the details matter!

> The one skill that separates bad programmers from good programmers is attention to detail. 
>
> Zed Shaw, _Learn Python the Hard Way_

In a written natural language, there are many ways to express the same idea. To 
make the consumption of information easier, people define style guides to enforce 
particularly effective ways of writing. This is even more true in coding; 
consistent style choices make scripts much easier to read. They become absolutely 
essential as projects become large (>1 person).

Some programming languages (\*cough\* *Java*) have multiple competing standards, 
and it's easy to imagine how messy this can get. Luckily, Python doesn't have 
this issue. The official standard, [PEP8](http://www.python.org/dev/peps/pep-0008/), 
is used everywhere. Unless you plan on hiding all the code you write from the 
outside world, you should learn it.

To help out coders, there are tools to test for compliance. The aptly named 
`pep8` library shows you where your code deviates from the PEP8 standard, and
`autopep8` goes a step further by trying to fix all of your errors for you.
These are both run from the shell, as

```
pep8 my_code.py
autopep8 my_code.py > my_new_code.py
```

These libraries won't always pick up everything, sadly. Furthermore, due to
the desire to maintain backward compatibility, there is some wiggle room in PEP8 
(see [this powerpoint](www.python.org/doc/essays/ppt/regrets/PythonRegrets.ppt) 
of Python regrets, made by the creator of the language). Here are some additional
rules to remember:

**PEP8 conventions missed by the `pep8` checker**

 * Variables and functions should be named  in `snake_case`. No capital letters.
   Classes are named in `UpperCamelCase`.
 * Multiline comments use `"""`, not `'''`.
 * Private methods and variables should be prefixed with an underscore, ie. 
   `_my_private_method()`
 
**Special rules outside of PEP8**

 * _Never_ use tabs. Ever.
 * Use list comprehensions over `map()`, `reduce()`, and `filter()`.
 * Avoid iterating through lists by index whenever possible.
 * Lambda functions should not be saved to a variable.
 
These rules might seem random (probably because they are), but, trust me: they
make collaborative coding _so much_ easier.

## Debuggers: for the deep-rooted errors

Linting will only catch the really obvious errors. For more complex issues,
(ie. bugs), you're going to want to follow the code's logic line by line. One
lazy way to do this is to put `print` statements everywhere, which allows you
to view variables over time. However, this gets messy quickly, and you lose
control of what variables you can see once you start executing.

This is where the Python DeBugger, or _pdb_, comes into play. With it, you can 
step through your code and watch as variables are changed.

All you have to do to use this is import the `pdb` module and call the 
`set_trace()` method.

```python
import pdb

# [ ... ]
# Your code here
# [ ... ]

pdb.set_trace()
```

Now, when you run the code, it will stop at whatever line you put `set_trace()`.
You'll be prompted to give a command. Some common commands include:

 * `continue` continues on to the next time a `set_trace()` line is hit
 * `print \*variable\*` prints the current value of a specified variable
 * `list` shows the source code around the `set_trace()` line
 * `args` prints the values of all the arguments in the current function

There are a lot more options, which can be found [here](http://docs.python.org/2/library/pdb.html), 
but these few should be enough to get you running with pdb.
 
## Profiling: making code fast

So, you've found your errors, those deep-rooted bugs, and even standardized
your code to conform to PEP8. But, for some reason, it's still really slow.
What can we do about this?

The first idea you might have is to time your code. Analogous to the `print`
statement debugging above, you could write some logic to print run times at
various points in your script.

```python
from time import time

t0 = time()
# run your code
print time() - t0
```

You can also time your entire script with:

```
time python my_code.py
```

While these both work, they're either too messy or not detailed enough. What we
really want is a breakdown of how long the computer spends running each part
of our code.

_Profilers_ provide a way to do just this. With Python, run your script in a
shell with this command

```
python -m cProfile -s time my_code.py
```

This returns the amount of functions called in the execution, along with a
breakdown of the time each function took. The `-s time` part sorts the output
by the time taken (which is usually what you care about). A sample output looks 
like

```
   2530004 function calls in 0.789 seconds
   Ordered by: internal time
   
   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
  1980000    0.253    0.000    0.253    0.000 my_code.py:23(<genexpr>)
    10000    0.190    0.000    0.780    0.000 my_code.py:16(evolve)
   220000    0.182    0.000    0.436    0.000 {sum}
   270000    0.150    0.000    0.150    0.000 my_code.py:28(neighbors)
        1    0.009    0.009    0.789    0.789 my_code.py:9(my_func)
    50000    0.004    0.000    0.004    0.000 {method 'add' of 'set' objects}
        1    0.000    0.000    0.000    0.000 {range}
```
With this information, you can go back into your code and adjust the logic to 
improve the speed bottlenecks.

## Segfaults: the scourge of C

Segmentation faults, abbreviated _segfaults_, are the worst debug errors in existence. 
Segfaults occur when the program tries to access a part of memory that it 
expects to be able to get to, and for whatever reason it is not available.
Because of this, segfaults only occur at runtime, so many tools won't even
catch it. To make it worse, a code with a segfault will return the most useless
error in the history of errors:

```
Segmentation fault
```

There's no traceback, and no explanation at all. Nothing.

Luckily, Python (generally) handles the memory issues that could cause a
segfault. Unluckily, many Python scripts interface with lower-level libraries;
these can segfault for, like, no reason.

If you happen upon a segfault in Python, you have two courses of action. If the
segfault is rare and affects a nonvital part of the code, then you can wrap it
in a standard Python exception and handle it in a regular manner. For this, the
[faulthandler](https://github.com/haypo/faulthandler/wiki/) library is useful.
If, however, you absolutely need the function, you're going to have to dive
into the C code and fix it yourself.

The tool you'll need to use is [valgrind](http://valgrind.org/), which is a
debugger + profiler + memory leak checker for compiled C code. It's so important
in C that tutorials in the language often include it within the first few
lessons. 

To use it, first compile the code of interest

```
g++ myCode.cc -o myCode
```

Then, run the compiled code through valgrind

```
valgrind ./myCode
```

Valgrind has a lot of features, each of which can be toggled in the command-line
call, but that's beyond the scope of this tutorial. Keep in mind that these
things exist, and, if you ever find yourself reading through some C, do yourself
a favor and go through a tutorial on all of this. I recommend [Learn C the Hard
Way](http://c.learncodethehardway.org/book/), but feel free to use whatever works.

## In conclusion

If you're new to coding, your head's probably already spinning with entirely
too much new information. That's okay. Remember that everything in this
tutorial is supplemental in nature; learning how to code is much more important
than knowing the tools that can help you code faster. However, keep in mind
that these tools all exist, and make it a goal to eventually come back here and
re-learn them with a clear mind.

If you're an experienced coder, these tools are exactly what you need to up
your game. Debuggers, profilers, and linters all save you valuable time, and
the PEP8 checker will really help in collaborative projects. Get used to them,
and the time investment will pay off.
