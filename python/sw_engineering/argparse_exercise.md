## Exercise 1

In the episode on [Randomness](http://software-carpentry.org/4_0/invperc/random)
we discussed how it is important (and saves time!) to use well-tested library
routines to generate random numbers.  This principle of reuse is also true for
any other task you need to carry out.  Python is packaged with an [extensive
library of modules](http://docs.python.org/library/index.html). What module
could you use to...
	
* parse a string containing a date or time (e.g. "Thursday, 27 May, 2010")? 
   > [`datetime`](http://docs.python.org/library/datetime.html) and in particular,
   > the function [`datetime.strptime()`](http://docs.python.org/library/datetime.html#datetime.datetime.strptime)


* inspect a bunch of files in a folder or subfolders?
   > One module for doing this is the
   > [`os.path`](http://docs.python.org/library/os.path.html).  It has lots of
   > functions for manipulating path names, and also this gem of a function:
   > [`os.path.walk`](http://docs.python.org/library/os.path.html#os.path.walk).

* manage command line arguments to your program?
   > [`argparse`](http://docs.python.org/library/argparse.html)

* access data on the web?
   > For very basic operations like downloading files, check out:
   > [`urllib`](http://docs.python.org/library/urllib.html)


## Exercise 2

In this lecture we wrote roughly the following snippet to parse the command line arguments.
    
```python
import sys

def fail(message):
 print message
 sys.exit(1)

def parse_arguments(arguments):
  '''Parse strings to get controlling parameters.'''

  try:
    grid_size   = int(arguments[0])
    value_range = int(arguments[1])
    random_seed = int(arguments[2])
  except IndexError:
    fail('Expected 3 arguments, got %d' % len(arguments))
  except ValueError:
    fail('Expected int arguments, got %s' % str(arguments))

  return grid_size, value_range, random_seed

if __name__ == '__main__':
 arguments = sys.argv[1:]
 grid_size, value_range, random_seed = parse_arguments(arguments)

  # print out the arguments
  print "grid size = %i \t value range = %i \t seed = %i" % \
     (grid_size, value_range, random_seed)
```

As we learned in Exercise 1, python comes with a library,
[`argparse`](http://docs.python.org/library/argparse.html), to do this in a way
that's easier to extend and with better error messages.  

In this exercise, you'll rework the above code to use the argparse library.
You'll need to import the `argparse` module, and rewrite `parse_arguments` to
use it.

**A hint on getting started:**

> You'll need to create an ArgumentParser object, and then call it's
> `add_argument` method for each of the arguments.


**A hint on setting up `argparse`:**

> The arguments are positional arguments so when you call `add_argument` you just
> need to supply a
> [name](http://docs.python.org/library/argparse.html#name-or-flags), and a type
> (extra hint: `int`).  Passing in a help message, default value, and so on is not
> mandatory.

**A hint if you're really stuck:**

> ```python
> parser.add_argument('grid_size', type=int, help="Grid size")
> ```

**A hint on accessing the parsed arguments:**

> If you created an argument named "foo", you can access it by calling:
> 
> ```python
> args = parser.parse_args()
> print args.foo
> ```

**Our answer:**
> ```python
> import sys, argparse
> 
> def parse_arguments(arguments):
>   '''Parse strings to get controlling parameters.'''
> 
>   parser = argparse.ArgumentParser()
>   parser.add_argument('grid_size',   type=int, help='Grid size')
>   parser.add_argument('value_range', type=int, help='Value range')
>   parser.add_argument('random_seed', type=int, help='Random seed')
>   args   = parser.parse_args(arguments)
> 
>   return args.grid_size, args.value_range, args.random_seed
> 
> if __name__ == '__main__':
>   arguments = sys.argv[1:]
>   grid_size, value_range, random_seed = parse_arguments(arguments)
>   print "grid size = %i \t value range = %i \t seed = %i" % \
>     (grid_size, value_range, random_seed)
> ```
>
> The `fail()` function is no longer necessary since `argparse.parse_args`
> exits and prints a help message if incorrect arguments are passed in.  Try
> running your program with only the argument "-h" to get a more verbose help
> message.
