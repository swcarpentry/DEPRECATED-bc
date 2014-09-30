---
layout: lesson
root: ../..
title: Simple Patterns
---

### Reading in the data
Let's start by reading data from two files named `notebook-1.txt` and `notebook-2.txt`. 
For each file, we will discard the header line, and keep the first 6 lines after it:

~~~
readings = []
for filename in ('notebook-1.txt', 'notebook-2.txt'):
    lines = open(filename, 'r').read().strip().split('\n')
    readings += lines[1:7] # We are ignoring the header line (lines[0]) here.

for r in readings:
    print r
~~~
{:class="in"}

This puts 6 lines from the first data file and 6 from the second
into the list `readings`:

~~~
Baker 1 2009-11-17      1223.0
Baker 1 2010-06-24      1122.7
Baker 2 2009-07-24      2819.0
Baker 2 2010-08-25      2971.6
Baker 1 2011-01-05      1410.0
Baker 2 2010-09-04      4671.6
Davison/May 23, 2010/1724.7
Pertwee/May 24, 2010/2103.8
Davison/June 19, 2010/1731.9
Davison/July 6, 2010/2010.7
Pertwee/Aug 4, 2010/1731.3
Pertwee/Sept 3, 2010/4981.0
~~~
{:class="out"}

Each element of `readings` is a *record* that the grad students have created. 
We will be testing our regular expressions against each record to see how well
we are matching the different record styles/formats as we go along.

### Pattern matching
As defined in the [Introduction](01-intro-regex.html), a regular expression is just a pattern that can match a string. In this section
we will aim to derive a regular expression which matches the records (strings of data) based on the date field. This will allow us to,
for example, list all the records that were created in May.

Without the use of regular expressions, we can determine if a record `r` contains the month
"06" using the `in` keyword (i.e. we are essentially asking the question 'does the record `r` contain the string "06"?'):

~~~
for r in readings:
    if '06' in r:
        print r
~~~
{:class="in"}
~~~
Baker 1 2010-06-24      1122.7
~~~
{:class="out"}

If we want to select all the records for two months we also have to use the `or` keyword:

~~~
for r in readings:
    if ('06' in r) or ('07' in r):
        print r
~~~
{:class="in"}
~~~
Baker 1 2010-06-24      1122.7
Baker 2 2009-07-24      2819.0
~~~
{:class="out"}

But if we say `'05' in r` it can match against the day "05" as well
as the month "05". This is not the desired behaviour we want. 
We can try to write a more complicated test that only
looks for the two-digit month in a particular place in the string, but
let's try using a regular expression to do this instead.

We will work up to our solution which uses a regular expression in stages.
We start by importing the regular expressions library, called `re`, then examine each record in the list `readings`. We will then
use the function `re.search` to try and find a match for the string `'06'` in a given record. 
If the function is successful in finding a match, then we print that record out.

~~~
import re
for r in readings:
    if re.search('06', r):
        print r
~~~
{:class="in"}
~~~
Baker 1 2010-06-24      1122.7
~~~
{:class="out"}

Note that the first argument to `re.search` is the pattern we are searching for,
written as a string. The second argument is the data we are searching
in. It's easy to reverse these accidentally, i.e., to put the data first
and the pattern second. This can be hard to track down, so please be
careful.

So far, our code that uses `re.search` does the same thing as `'06' in r`. But if we want to match
`'06'` or `'07'`, regular expressions let us combine the two comparisons
in a single expression:

~~~
for r in readings:
    if re.search('06|07', r):
        print r
~~~
{:class="in"}
~~~
Baker 1 2010-06-24      1122.7
Baker 2 2009-07-24      2819.0
~~~
{:class="out"}

The vertical bar in the pattern means "or". It tells regular expression
library that we want to match either the specified text on the left of the vertical bar, *or* the text
on the right of the vertical bar.

We are going to be throwing a lot of regular expressions against our
data, so to help us test whether different regular expressions are doing the right thing,
let's write a function that will tell us which records match a
particular pattern. Our function `show_matches` takes a pattern and a
list of strings as arguments. It prints out two stars as a marker if the
pattern matches a string, and just indents with blanks if it does not:

~~~
def show_matches(pattern, strings):
    for s in strings:
        if re.search(pattern, s):
            print '**', s
        else:
            print '  ', s
~~~
{:class="in"}

If we use this function to match `'06|07'` against the data we read in
earlier, it prints stars beside the two records that have month `'06'`
or month `'07'`:

~~~
show_matches('06|07', readings)
~~~
{:class="in"}
~~~
   Baker 1  2009-11-17  1223.0
** Baker 1  2010-06-24  1122.7
** Baker 2  2009-07-24  2819.0
   Baker 2  2010-08-25  2971.6
   Baker 1  2011-01-05  1410.0
   Baker 2  2010-09-04  4671.6
   Davison/May 23, 2010/1724.7
   Pertwee/May 24, 2010/2103.8
   Davison/June 19, 2010/1731.9
   Davison/July 6, 2010/2010.7
   Pertwee/Aug 4, 2010/1731.3
   Pertwee/Sept 3, 2010/4981.0
~~~
{:class="out"}

But if we change the pattern `'06|7'` (without a '0' in front of the
'7'), the pattern seems to match a lot of things that don't have the
month `'06'` or `'07'`:

~~~
show_matches('06|7', readings)
~~~
{:class="in"}
~~~
** Baker 1  2009-11-17  1223.0
** Baker 1  2010-06-24  1122.7
** Baker 2  2009-07-24  2819.0
** Baker 2  2010-08-25  2971.6
   Baker 1  2011-01-05  1410.0
** Baker 2  2010-09-04  4671.6
** Davison/May 23, 2010/1724.7
   Pertwee/May 24, 2010/2103.8
** Davison/June 19, 2010/1731.9
** Davison/July 6, 2010/2010.7
** Pertwee/Aug 4, 2010/1731.3
   Pertwee/Sept 3, 2010/4981.0
~~~
{:class="out"}

To understand why, think back to mathematics. The expression *ab+c*
means "a times b, plus c" because multiplication has higher precedence
than addition. If we want to force the other meaning, we have to use
parentheses and write *a(b+c)*.

The same is true for regular expressions. Adjacency has higher
precedence than "or", so the pattern `'06|7'` means, "Either `'06'` or
the digit `'7'`". If we look back at our data, there are a lot of 7's in
our file, and this pattern is matching all of them.

If we want to match `'06'` or `'07'` without repeating the digit '0', we
have to parenthesize it as `'0(6|7)'`. Having said that, most people
probably find the expression `'06|07'` more readable anyway. Note that
the expression inside the parentheses is a regular expression in its own right,
and is therefore referred to as a *sub-expression*.

Let's go back to our function and our data. If we use the pattern
`'05'`, then as we said earlier, we will match records that have '05' as
the day as well as those with '05' as the month. We can force our match
to do the right thing by taking advantage of context. If the date is
formatted as YYYY-MM-DD then there should be a dash `'-'` before and
after the month, but only before the day. The pattern `'-05-'` should
therefore only match a month of '05'. Sure enough, if we give this
pattern to our function it doesn't match any records. This is the
correct answer, since we don't have any readings in this sample of our
data set for May.

### Extracting data
Matching is useful, but what we really want to do is extract the year,
the month, and the day from our data so that we can reformat them.
Parentheses can help here too: when a regular expression matches a piece
of text, the library automatically remembers what matched against every
parenthesized sub-expression.

Here's a simple example:

~~~
match = re.search('(2009|2010|2011)',
                   'Baker 1\t2009-11-17\t1223.0')
print match.group(1)
~~~
{:class="in"}

The first string is our pattern. It will match 2009, 2010, or 2011, and
the parentheses around it will make the library remember which of those
three strings was matched. The second string is just the first record
from our data. (Remember, `'\t'` represents a tab.)

When `re.search` is called, it returns `None` if it doesn't find a
match, or a special [match object](https://docs.python.org/2/library/re.html#match-objects) if it
does. The call to `match.group` returns the text that matched the
sub-expression inside the specified set of parentheses, counting from the
left. Since this pattern only has one set of parentheses,
`match.group(1)` returns whatever text matched what's inside them.

The way sub-expressions are numbered sometimes trips people up. While
Python normally counts from 0, the first match in a regular expression
is extracted with `match.group(1)`, the second with 2, and so forth. The
reason is that `match.group(0)` returns all of the text that the entire
pattern matched.

### The 'dot' character
What if we want to match the month as well as the year? A regular
expression to match legal months would be
`'(01|02|03|04|05|06|07|08|09|10|11|12)'`. An expression to match days
would be three times longer, which would be hard to type and (more
importantly) hard to read.

Instead, we can use the dot character `'.'` to match any single
character. For example, the expression `'....'` matches exactly
four characters, and `'....-..-..'` matches four characters, a dash, two
more characters, another dash, and two more characters. If we put each
set of dots in parentheses as `'(....)-(..)-(..)'` the three groups
should record the year, month, and day each time there's a successful
match.

Let's test that out by calling `re.search` with the pattern we just
described and the first record from our data:

~~~
match = re.search('(....)-(..)-(..)',
                   'Baker 1\t2009-11-17\t1223.0')
print match.group(1), match.group(2), match.group(3)
~~~
{:class="in"}
~~~
2009 11 17
~~~
{:class="out"}

When we print out the three groups, we get `'2009'`, `'11'`, and `'17'`,
just as we wanted. Try doing *that* with substring searches...

**Tip:** If you want to match an actual dot/period/full-stop '.' character, you must place a single backward-slash before it (i.e. `'\.'`). 
Just using a single dot on its own (i.e. using `'.'`) will result in *any* single character being matched, as demonstrated in the section above.

<div class="keypoints" markdown="1">

#### Key Points
* A regular expression is just a pattern that can match a string. 
* A single alphanumeric character (or string of alphanumeric characters like `'05'` or `'06'`) is a regular expression in its own right, and only matches itself. For example, `'06'` will not match the string "6", nor will it match "60"; it will only match "06". Also, `'A'` only matches an upper-case A, but not a lower-case one. 
* The vertical bar character `'|'` means "or", and can be used to combine two comparisons in just one regular expression. 
* A dot `'.'` character matches any single character, and we use parentheses to enforce grouping and to remember things.
* Build a pattern by starting with something simple that matches part of the data we're working with, then add to it piece by piece. 
* Test the pattern against our data each time we make a change, but also test that it *doesn't* match things that it shouldn't, because false positive can be very hard to track down.

</div>

<div class="challenge" markdown="1">
Write a regular expression that matches each of the following binary strings: `000`, `101`, and `111`.
</div>

<div class="challenge" markdown="1">
Write a regular expression to match all binary strings that are at least 1 digit long and at most 4 digits long.
</div>

<div class="challenge" markdown="1">
Write a program which reads in a file containing the following words: `hello, working, telling, as, meaningful, cold, world, caring, ingrid`, and uses a regular expression to match all the words that end in 'ing'.
</div>
