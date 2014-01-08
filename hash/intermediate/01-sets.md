---
layout: lesson
root: ../..
title: Sets
level: intermediate
---
<div class="objectives" markdown="1">
## Objectives

*   Explain why some programs that use lists become proportionally
    slower as data sizes increase.
*   Explain the three adjectives in "unordered collection of distinct
    values".
*   Use a set to eliminate duplicate values from data.
</div>

## Lesson

Let's start with something simpler than our actual inventory problem.
Suppose we have a list of all the atoms in the warehouse, and we want to
know which different kinds we have—not how many, but just their types.
We could solve this problem using a list to store the unique atomic
symbols we have seen. Here's a function to add a new atom to the list:

    def another_atom(seen, atom):
        for i in range(len(seen)):
            if seen[i] == atom:
                return # atom is already present, so do not re-add
        seen.append(atom)

`another_atom`'s arguments are a list of the unique atoms we've already
seen, and the symbol of the atom we're adding. Inside the function, we
loop over the atoms that are already in the list. If we find the one
we're trying to add, we exit the function immediately: we aren't
supposed to have duplicates in our list, so there's nothing to add. If
we reach the end of the list without finding this symbol, though, we
append it. This is a common [design
pattern](glossary.html#design-pattern): either we find pre-existing data
in a loop and return right away, or take some default action if we
finish the loop without finding a match.

Let's watch this function in action. We start with an empty list. If the
first atomic symbol is `'Na'`, we find no match (since the list is
empty), so we add it. The next symbol is `'Fe'`; it doesn't match
`'Na'`, so we add it as well. Our third symbol is `'Na'` again. It
matches the first entry in the list, so we exit the function
immediately.

  Before           Adding   After
  ---------------- -------- ----------------
  `[]`             `'Na'`   `['Na']`
  `['Na']`         `'Fe'`   `['Na', 'Fe']`
  `['Na', 'Fe']`   `'Na'`   `['Na', 'Fe']`

This code works, but it is inefficient. Suppose there are *V* distinct
atomic symbols in our data, and *N* symbols in total. Each time we add
an observation to our list, we have to look through an average of *V/2*
entries. The total running time for our program is therefore
approximately *NV/2*. If *V* is small, this is only a few times larger
than *N*, but what happens if we're keeping track of something like
patient records rather than atoms? In that case, most values are
distinct, so *V* is approximately the same as *N*, which means that our
running time is proportional to *N^2^/2*. That's bad news: if we double
the size of our data set, our program runs four times slower, and if we
double it again, our program will have slowed down by a factor of 16.

There's a better way to solve this problem that is simpler to use and
runs much faster. The trick is to use a [set](glossary.html#set) to
store the symbols. A set is an unordered collection of distinct items.
The word "collection" means that a set can hold zero or more values. The
word "distinct" means that any particular value is either in the set or
not: a set can't store two or more copies of the same thing. And
finally, "unordered" means that values are simply "in" the set. They're
not in any particular order, and there's no first value or last value.
(They actually are stored in some order, but as we'll discuss in [the
next section](#s:storage), that order is as random as the computer can
make it.)

To create a set, we simply write down its elements inside curly braces:

    >>> primes = {3, 5, 7}

![A Simple Set](setdict/simple_set.png)

Figure 1: A Simple Set

However, we have to use `set()` to create an empty set, because the
symbol `{}` was already being used for something else when sets were
added to Python:

    >>> even_primes = set() # not '{}' as in math

We'll meet that "something else" [later in this chapter](#s:dict).

To see what we can do with sets, let's create three holding the integers
0 through 9, the first half of that same range of numbers (0 through 4),
and the odd values 1, 3, 5, 7, and 9:

    >>> ten  = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
    >>> lows = {0, 1, 2, 3, 4}
    >>> odds = {1, 3, 5, 7, 9}

If we ask Python to display one of our sets, it shows us this:

    >>> print lows
    set([0, 1, 2, 3, 4])

rather than using the curly-bracket notation. I personally regard this
as a design flaw, but it does remind us that we can create always create
a set from a list.

Sets have methods just like strings and lists, and, like the methods of
strings and lists, most of those methods create new sets instead of
modifying the set they are called for. These three come straight from
mathematics:

    >>> print lows.union(odds)
    set([0, 1, 2, 3, 4, 5, 7, 9])
    >>> print lows.intersection(odds)
    set([1, 3])
    >>> print lows.difference(odds)
    set([0, 2, 4])

Another method that creates a new set is `symmetric_difference`, which
is sometimes called "exclusive or":

    >>> print lows.symmetric_difference(odds)
    set([0, 2, 4, 5, 7, 9])

It returns the values that are in one set or another, but not in both.

Not all set methods return new sets. For example, `issubset` returns
`True` or `False` depending on whether all the elements in one set are
present in another:

    >>> print lows.issubset(ten)
    True

A complementary method called `issuperset` also exists, and does the
obvious thing:

    >>> print lows.issuperset(odds)
    False

We can count how many values are in a set using `len` (just as we would
to find the length of a list or string), and check whether a particular
value is in the set or not using `in`:

    >>> print len(odds)
    7
    >>> print 6 in odds
    False

Some methods modify the sets they are called for. The most commonly used
is `add`, which adds an element to the set:

    >>> lows.add(9)
    >>> print lows
    set([0, 1, 2, 3, 4, 9])

If the thing being added is already in the set, `add` has no effect,
because any specific thing can appear in a set at most once:

    >>> lows.add(9)
    >>> print lows
    set([0, 1, 2, 3, 4, 9])

This behavior is different from that of `list.append`, which always adds
a new element to a list.

Finally, we can remove individual elements from the set:

    >>> lows.remove(0)
    >>> print lows
    set([1, 2, 3, 4])

or clear it entirely:

    >>> lows.clear()
    >>> print lows
    set()

Removing elements is similar to deleting things from a list, but there's
an important difference. When we delete something from a list, we
specify its *location*. When we delete something from a set, though, we
must specify the *value* that we want to take out, because sets are not
ordered. If that value isn't in the set, `remove` does nothing.

To help make programs easier to type and read, most of the methods we've
just seen can be written using arithmetic operators as well. For
example, instead of `lows.issubset(ten)`, we can write `lows <= ten`,
just as if we were using pen and paper. There are even a couple of
operators, like the strict subset test `<`, that don't have long-winded
equivalents.

  Operation           As Method                           Using Operator
  ------------------- ----------------------------------- ----------------
  *difference*        `lows.difference(odds)`             `lows - odds`
  *intersection*      `lows.intersection(odds)`           `lows & odds`
  *subset*            `lows.issubset(ten)`                `lows <= ten`
  *strict subset*                                         `lows < ten`
  *superset*          `lows.issuperset(ten)`              `lows >= odds`
  *strict superset*                                       `lows >= odds`
  *exclusive or*      `lows.symmetric_difference(odds)`   `lows ^ odds`
  *union*             `lows.union(odds)`                  `lows | odds`

The fact that the values in a set are distinct makes them a convenient
way to get rid of duplicate values, like the "unique atoms" problem at
the start of this section. Suppose we have a file containing the names
of all the atoms in our warehouse, and our task is to produce a list of
the their types. Here's how simple that code is:

~~~~ {src="setdict/unique_atoms.py"}
import sys
filename = sys.argv[1]
source = open(filename, 'r')
atoms = set()
for line in source:
    name = line.strip()
    atoms.add(name)
source.close()
print atoms
~~~~

We start by opening the file and creating an empty set which we will
fill with atomic symbols. As we read the lines in the file, we strip off
any whitespace (such as the newline character at the end of the line)
and put the resulting strings in the set. When we're done, we print the
set. If our input is the file:

        Na
        Fe
        Na
        Si
        Pd
        Na

then our output is:

    set(['Fe', 'Si', 'Na'])

The right atoms are there, but what are those extra square brackets for?
The answer is that if we want to construct a set with values using
`set()`, we have to pass those values in a single object, such as a
list. This syntax:

    set('Na', 'Fe', 'Fl')

does *not* work, even though it seems more natural. On the other hand,
this means that we can construct a set from almost anything that a `for`
loop can iterate over:

    >>> set('lithium')
    set(['i', 'h', 'm', 'l', 'u', 't'])

But hang on: if we're adding characters to the set in the order `'l'`,
`'i'`, `'t'`, `'h'`, `'i'`, `'u'`, `'m'`, why does Python show them in
the order `'i'`, `'h'`, `'m'`, `'l'`, `'u'`, `'t'`? To answer that
question, we need to look at how sets are actually stored, and why
they're stored that way.

<div class="keypoints" markdown="1">
## Key Points

*   Use sets to store distinct unique values.
*   Create sets using `set()` or `{v1, v2, ...}`.
*   Sets are mutable, i.e., they can be updated in place like lists.
*   A loop over a set produces each element once, in arbitrary order.
*   Use sets to find unique things.
</div>

<div class="challenges" markdown="1">
## Challenges

1.  Mathematicians are quite comfortable negating sets: for example, the
    negation of the set `{1, 2}` is all numbers that aren't 1 or 2. Why
    don't Python's sets have a `not` operator?
2.  Fan has created a set containing the names of five noble gases:

        >>> print gases
        set(['helium', 'argon', 'neon', 'xenon', 'radon'])

    He would like to print them in alphabetical order. What is one
    simple way to do this? (Hint: the `list` function converts its
    arguments to a list.)

3.  Fan has the following code:

        left = {'He', 'Ar', 'Ne'}
        right = set()
        while len(left) > len(right):
            temp = left.pop()
            right.add(temp)

    What values could \`left\` and \`right\` have after this code is
    finished running? Explain why your answer makes this code hard to
    test.

4.  Fan has written the following code:

        left = {'He', 'Ar', 'Ne'}
        right = {'Ar', 'Xe'}
        for element in left:                # X
            if element not in right:        # X
                right.add(element)          # X
        assert left.issubset(right)

    What single line could be used in place of the three marked with 'X'
    to achieve the same effect?

5.  Fan has written a program to print the names of the distinct atoms
    in a data file:

    ~~~~ {src="setdict/print_names.py"}
    # Print the name of each atom in the data file once.
    reader = open('atoms.txt', 'r')
    seen = set()
    for line in reader:
        name = line.strip()
        if name in seen:
            print name
        else:
            seen.add(name)
    reader.close()
    ~~~~

    When he runs the program on this data file:

    ~~~~ {src="setdict/atoms.txt"}
    Na
    Fe
    Na
    ~~~~

    it only prints:

        Na

    What is the simplest change you can make to the program so that it
    produces the correct answer?
</div>
