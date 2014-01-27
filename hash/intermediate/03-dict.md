---
layout: lesson
root: ../..
title: Dictionaries
level: intermediate
---
<div class="objectives" markdown="1">
## Objectives

*   Explain the similarities and differences between sets and
    dictionaries.
*   Perform common operations on dictionaries.
</div>

## Lesson

Now that we know how to find out what kinds of atoms are in our
inventory, we want to find out how many of each we have. Our input is a
list of several thousand atomic symbols, and the output we want is a
list of names and counts.

Once again, we could use a list to store names and counts, but the right
solution is to use another new data strucure called a
[dictionary](glossary.html#dictionary). A dictionary is a unordered
collection of key-value pairs ([Fixture 5](#f:simple_dict)). The keys
are immutable, unique, and unordered, just like the elements of a set.
There are no restrictions on the values stored with those keys: they
don't have to be immutable or unique. However, we can only look up
entries by their keys, not by their values.

![A Simple Dictionary](setdict/simple_dict.png)

Figure 5: A Simple Dictionary

We create a new dictionary by putting key-value pairs inside curly
braces with a colon between the two parts of each pair:

    >>> birthdays = {'Newton' : 1642, 'Darwin' : 1809}

The dictionary's keys are the strings `'Newton'` and `'Darwin'`. The
value associated with `'Newton'` is 1642, while the value associated
with `'Darwin'` is 1809. We can think of this as a two-column table:

  Key          Value
  ------------ -------
  `'Newton'`   1642
  `'Darwin'`   1809

but it's important to remember that the entries aren't necessarily
stored in this order (or any other specific order).

We can get the value associated with a key by putting the key in square
brackets:

    >>> print birthdays['Newton']
    1642

This looks just like subscripting a string or list, except dictionary
keys don't have to be integers—they can be strings, tuples, or any other
immutable object. It's just like using a phonebook or a real dictionary:
instead of looking things up by location using an integer index, we look
things up by name.

If we want to add another entry to a dictionary, we just assign a value
to the key, just as we create a new variable in a program by assigning
it a value:

    >>> birthdays['Turing'] = 1612
    >>> print birthdays
    {'Turing' : 1612, 'Newton' : 1642, 'Darwin' : 1809}

If the key is already in the dictionary, assignment replaces the value
associated with it rather than adding another entry (since each key can
appear at most once). Let's fix Turing's birthday by replacing 1612 with
1912:

    >>> birthdays['Turing'] = 1912
    >>> print birthdays
    {'Turing' : 1912, 'Newton' : 1642, 'Darwin' : 1809}

Trying to get the value associated with a key that *isn't* in the
dictionary is an error, just like trying to access a nonexistent
variable or get an out-of-bounds element from a list. For example, let's
try to look up Florence Nightingale's birthday:

    >>> print birthdays['Nightingale']
    KeyError: 'Nightingale'

If we're not sure whether a key is in a dictionary or not, we can test
for it using `in`:

    >>> print 'Nightingale' in birthdays
    False
    >>> print 'Darwin' in birthdays
    True

And we can see how many entries are in the dictionary using `len`:

    >>> print len(birthdays)
    3

and loop over the keys in a dictionary using `for`:

    >>> for name in birthdays:
    ...     print name, birthdays[name]
    ...
    Turing 1912
    Newton 1642
    Darwin 1809

This is a little bit different from looping over a list. When we loop
over a list we get the values in the list. When we loop over a
dictionary, on the other hand, the loop gives us the keys, which we can
use to look up the values.

We're now ready to count atoms. The main body of our program looks like
this:

~~~~ {src="setdict/count_atoms.py"}
import sys

if __name__ == '__main__':
    reader = open(sys.argv[1], 'r')
    lines = reader.readlines()
    reader.close()
    count = count_atoms(lines)
    for atom in count:
        print atom, count[atom]
~~~~

The first three lines read the input file into a list of strings. We
then call a function `count_atoms` to turn that list into a dictionary
of atomic symbols and counts. Once we have that dictionary, we use a
loop like the one we just saw to print out its contents.

Here's the function that does the counting:

~~~~ {src="setdict/count_atoms.py"}
def count_atoms(lines):
  '''Count unique atoms, returning a dictionary.'''

    result = {}
    for atom in lines:
        atom = atom.strip()
        if atom not in result:
            result[atom] = 1
        else:
            result[atom] = result[atom] + 1

    return result
~~~~

We start with a docstring to explain the function's purpose to whoever
has to read it next. We then create an empty dictionary to fill with
data, and use a loop to process the lines from the input file one by
one. Notice that the empty dictionary is written `{}`: this is the
"[previous use](#a:previous-use)" we referred to when explaining why an
empty set had to be written `set()`.

After stripping whitespace off the atom's symbol, we check to see if
we've seen it before. If we haven't, we set its count to 1, because
we've now seen that atom one time. If we *have* seen it before, we add
one to the previous count and store that new value back in the
dictionary. When the loop is done, we return the dictionary we have
created.

Let's watch this function in action. Before we read any data, our
dictionary is empty. After we see `'Na'` for the first time, our
dictionary has one entry: its key is `'Na'`, and its value is 1. When we
see `'Fe'`, we add another entry to the dictionary with that string as a
key and 1 as a value. Finally, when we see `'Na'` for the second time,
we add one to its count.

  Input     Dictionary
  --------- ------------------------
  *start*   `{}`
  `Na`      `{'Na' : 1}`
  `Fe`      `{'Na' : 1, 'Fe' : 1}`
  `Na`      `{'Na' : 2, 'Fe' : 1}`

Just as we use tuples for multi-part entries in sets, we can use them
for multi-part keys in dictionaries. For example, if we want to store
the years in which scientists were born using their full names, we could
do this:

    birthdays = {
        ('Isaac', 'Newton') : 1642,
        ('Charles', 'Robert', 'Darwin') : 1809,
        ('Alan', 'Mathison', 'Turing') : 1912
    }

If we do this, though, we always have to look things up by the full key:
there is no way to ask for all the entries whose keys contain the word
`'Darwin'`, because Python cannot match part of a tuple.

If we think of a dictionary as a two-column table, it is occasionally
useful to get one or the other column, i.e., just the keys or just the
values:

    all_keys = birthdays.keys()
    print all_keys
    [('Isaac', 'Newton'), ('Alan', 'Mathison', 'Turing'), ('Charles', 'Robert', 'Darwin')]
    all_values = birthday.values()
    print all_values
    [1642, 1912, 1809]

These methods should be used sparingly: the dictionary doesn't store the
keys or values in a list, these methods both actually create a new list
as their result. In particular, we *shouldn't* loop over a dictionary's
entries like this:

    for key in some_dict.keys():
        ...do something with key and some_dict[key]

since "`for key in some_dict`" is shorter and much more efficient.

### Summary

-   Use dictionaries to store key-value pairs with distinct keys.
-   Create dictionaries using `{k1:v1, k2:v2, ...}`
-   Dictionaries are mutable, i.e., they can be updated in place.
-   Dictionary keys must be immutable, but values can be anything.
-   Use tuples to store multi-part keys in dictionaries.
-   `dict[key]` refers to the dictionary entry with a particular key.
-   `key in dict` tests whether a key is in a dictionary.
-   `len(dict)` returns the number of entries in a dictionary.
-   A loop over a dictionary produces each key once, in arbitrary order.
-   `dict.keys()` creates a list of the keys in a dictionary.
-   `dict.values()` creates a list of the keys in a dictionary.

### Challenges

1.  What is one possible output of the following program? And why does
    this question say "one possible output" instead of "*the* output"?

        periods = {'Mercury' : 87.97, 'Venus' : 224.70}
        print periods
        periods.update({'Earth' : 3.6526, 'Mars' : 686.98})
        print periods
        periods['Earthy'] = 365.26
        print periods

2.  Fan has a table with the pH levels of samples as the keys, and the
    percentage of carbon-12 as the values:

      pH     C12
      ------ ------
      7.43   0.48
      7.51   0.47
      7.56   0.45

    He needs to interpolate between these values, i.e., to predict the
    percentage of carbon-12 in the sample for a pH 7.50. Will storing
    his data in a dictionary:

        {7.43 : 0.48, 7.51 : 0.47, 7.56 : 0.45}

    be any more efficient than storing it in a list of pairs:

        [ [7.43, 0.48], [7.51, 0.47], [7.56, 0.45] ]

    Why or why not?

3.  Before sets were added to Python, people frequently imitated them
    using dictionaries; the dictionary's keys were the set's elements,
    and the dictionary's values were all `None`. This function
    calculates the intersection of two such "sets":

        def setdict_intersect(left, right):
            '''Return new dictionary with intersection of keys from left and right.'''
            result = {}
            for key in left:
                if key in right:
                    result[key] = None
            return result

    Write a function `setdict_union` that calculates the union of two
    sets represented in this way.

4.  What does the following function do? Explain when and why you would
    use it, and write a small example that calls it with sample data.

        def show(writer, format, data):
            keys = data.keys()
            keys.sort()
            for k in keys:
                print >> writer, format % (key, data[key])

5.  Dictionaries are more general than lists, since you can trivially
    simulate a list like `['first', 'second', 'third']` using
    `{0 : 'first', 1 : 'second', 2 : 'third'}`. Given that, when and why
    should you use a list rather than a dictionary?


Summing Up
----------

Every programmer meets lists (or arrays or matrices) early in her
career. Many in science never meet sets and dictionaries, and that's a
shame: they often make programs easier to write and faster to run at the
same time.

Before we leave this topic, try running the function `globals` at an
interactive Python prompt:

    >>> globals()
    {'__builtins__': <module '__builtin__' (built-in)>,
     '__doc__': None,
     '__name__': '__main__',
     '__package__': None}

That's right—Python actually stores the program's variables in a
dictionary. In fact, it uses one dictionary for the global variables and
one for each currently-active function call:

    >>> def example(first, second):
    ...     print 'globals in example', globals()
    ...     print 'locals in example', locals()
    ... 
    >>> example(22, 33)
    globals in example {'__builtins__': <module '__builtin__' (built-in)>,
                        '__doc__': None,
                        '__name__': '__main__',
                        '__package__': None,
                        'example': <function example at 0x50b630>}
    locals in example {'second': 33,
                       'first': 22}

You now know everything you need to know in order to build a programming
language of your own. But please don't: the world will be much better
off if you keep doing science instead.

<div class="keypoints" markdown="1">
## Key Points

*   FIXME
</div>

<div class="challenges" markdown="1">
## Challenges

1.  FIXME
</div>
