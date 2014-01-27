---
layout: lesson
root: ../..
title: Storage
level: intermediate
---
<div class="objectives" markdown="1">
## Objectives

*   Draw a diagram showing how hash tables are implemented, and
    correctly label the main parts.
*   Explain the purpose of a hash function.
*   Explain why using mutable values as keys in a hash table can cause
    problems.
*   Correctly identify the error messages Python produces when programs
    try to put mutable values in hashed data structures.
*   Explain the similarities and differences between tuples and lists.
*   Explain why using tuples is better than concatenating values when
    storing multi-part data in hashed data structures.
</div>

## Lesson

Let's create a set, add the string `'lithium'` to it (as a single item,
not character by character), and print the result:

    >>> things = set()
    >>> things.add('lithium')
    >>> print things
    set(['lithium'])

As expected, the string is in the set. Now let's try adding a list to
the same set:

    >>> things.add([1, 2, 3])
    TypeError: unhashable type: 'list'

Why doesn't that work? And what does that word "unhashable" mean?

When we create a set, the computer allocates a block of memory to store
references to the set's elements. When we add something to the set, or
try to look something up, the computer uses a [hash
function](glossary.html#hash-function) to figure out where to look. A
hash function is any function that produces a seemingly-random number
when given some data as input. For example, one way to hash a string is
to add up the numerical values of its characters. If the string is
"zebra", those values are 97 for lower-case 'a', 98 for lower-case 'b',
and so on up to 122 for lower-case 'z'. When we add them up, we will
always get the same result: in this case, 532. If our hash table has 8
slots, we can take the remainder `532%8=4` to figure out where to store
a reference to our string in the hash table ([Figure
2](#f:set_storage_string)).

![Hashing a String](setdict/set_storage_string.png)

Figure 2: Hashing a String

Now let's take a look at how a list would be stored. If the list
contains the same five characters, so that its hash code is still 4, it
would be stored as shown in [Figure 3](#f:set_storage_list):

![Hashing a List](setdict/set_storage_list.png)

Figure 3: Hashing a List

But what happens if we change the characters in the list after we've
added it to the set? For example, suppose that we change the first
letter in the list from 'z' to 'X'. The hash function's value is now 498
instead of 532, which means that the modified list belongs in slot 2
rather than slot 4. However, the reference to the list is still in the
old location: the set doesn't know that the list's contents have
changed, so it hasn't moved its reference to the right location ([Figure
4](#f:set_storage_mutate)):

![After Mutation](setdict/set_storage_mutate.png)

Figure 4: After Mutation

This is bad news. If we now ask, "Is the list containing 's', 'e', 'b',
'r', and 'a' in the set?" the answer will be "no", because the reference
to the list isn't stored in the location that our hash function tells us
to look. It's as if someone changed their name from "Tom Riddle" to
"Lord Voldemort", but we left all the personnel records filed under 'R'.

This problem arises with any [mutable](glossary.html#mutable)
structure—i.e., any structure whose contents or value can be changed
after its creation. Integers and strings are safe to hash because their
values are fixed, but the whole point of lists is that we can grow them,
shrink them, and overwrite their contents.

Different languages and libraries handle this problem in different ways.
One option is to have each list keep track of the sets that it is in,
and move itself whenever its values change. However, this is expensive:
every time a program touched a list, it would have to see if it was in
any sets, and if it was, recalculate its hash code and update all the
references to it.

A second option is to shrug and say, "It's the programmer's fault." This
is what most languages do, but it's also expensive: programmers can
spend hours tracking down the bugs that arise from data being in the
wrong place.

Python uses a third option: it only allows programmers to put
[immutable](glossary.html#immutable) values in sets. After all, if
something's value can't change, neither can its hash code or its
location in a hash table.

But if sets can only hold immutable values, what do we do with mutable
ones? In particular, how should we store things like (x,y) coordinates,
which are naturally represented as lists, or people's names, which are
naturally represented as lists of first, middle, and last names? Again,
there are several options.

The first is to concatenate those values somehow. For example, if we
want to store "Charles" and "Darwin", we'd create the string "Charles
Darwin" and store that. This is simple to do, but our code will wind up
being littered with string joins and string splits, which will make it
slower to run and harder to read. More importantly, it's only safe to do
if we can find a concatenator that can never come up in our data. (If we
join "Paul Antoine" and "St. Cyr" using a space, there would be three
possible ways to split it apart again.)

The second option—the right one—is to use [tuples](glossary.html#tuple)
instead of lists. A tuple is an immutable list, i.e., a sequence of
values that cannot be changed after its creation. Tuples are created
exactly like lists, except we use parentheses instead of square
brackets:

    >>> full_name = ('Charles', 'Darwin')

They are indexed the same way, too, and functions like `len` do exactly
what we'd expect:

    >>> print full_name[0]
    Charles
    >>> len(full_name)
    2

What we *cannot* do is assign a new value to a tuple element, i.e.,
change the tuple after it has been created:

    >>> full_name[0] = 'Erasmus'
    TypeError: 'tuple' object does not support item assignment

This means that a tuple's hash code never changes, and *that* means that
tuples can be put in sets:

    >>> names = set()
    >>> names.add(('Charles', 'Darwin'))
    >>> print names
    set([('Charles', 'Darwin')])

<div class="keypoints" markdown="1">
## Key Points

*   Sets are stored in hash tables, which guarantee fast access for
    arbitrary values.
*   The values in sets must be immutable to prevent hash tables
    misplacing them.
*   Use tuples to store multi-part values in sets.
</div>

<div class="challenges" markdown="1">
## Challenges

1.  A friend of yours argues, "Finding a value in an unordered list of
    length *N* takes *N/2* steps on average. Finding it in a hash table
    takes only one step, but it's a more expensive step, since we have
    to calculate a hash code for that value. We should therefore use
    lists for small data sets, and only use things like sets for large
    ones." Explain the flaws in your friend's reasoning.
2.  Nelle has inherited the following function:

        def is_sample_repeated(left_channel, right_channel, history):
            '''Report repeated samples.  Both channels' values are integers in [0..10] inclusive.'''
            combined = 1000 * left_channel + right_channel
            if combined in history:
                return True
            else:
                history.add(combined)
                return False

    How would you improve this function, and why?

3.  Nelle has a function that extracts the latitudes and longitudes of
    data collection sites from a file:

        >>> sites = extract_sites('north-pacific.dat')
        >>> print sites[:3]
        [[52.097, -173.505], [52.071, -173.510], [51.985, -173.507]]

    Write another function called `filter_duplicate_sites` that takes a
    list of this kind as its only input, and returns a set (not a list)
    containing only the unique latitude/longitude pairs.

4.  A list containing just the number 5 is written as `[5]`, and a set
    containing just that same number is written as `{5}`. However, a
    tuple containing just that number must be written with a comma as
    `(5,)`. Why?
</div>
