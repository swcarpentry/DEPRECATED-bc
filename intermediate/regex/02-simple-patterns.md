# Simple Patterns


Let's start by reading data from two files, discarding the headers, and
keeping the first few lines of each:

    readings = []
    for filename in ('data-1.txt', 'data-2.txt'):
      lines = open(filename, 'r').read().strip().split('\n')
      readings += lines[2:8]

    for r in readings:
      print r

This puts six lines from the first data file and six from the second
into the list `readings`:

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

We will test our regular expressions against this data to see how well
we are matching different record formats as we go along.

Without regular expressions, we can select records that have the month
"06" using `if '06' in record`:

    for r in readings:
      if '06' in r:
        print r
    Baker 1 2010-06-24      1122.7

If we want to select data for two months we have to use
`if ('06' in record) or ('07' in record)`.

    for r in readings:
      if ('06' in r) or ('07' in r):
        print r
    Baker 1 2010-06-24      1122.7
    Baker 2 2009-07-24      2819.0

But if we say `'05' in record` it can match against the day "05" as well
as the month "05". We can try to write a more complicated test that only
looks for the two-digit month in a particular place in the string, but
let's try using a regular expression instead.

We will work up to our solution in stages. We start by importing the
regular expressions library, then examine each record in `readings`. If
a regular expression search can find a match for the string `'06'` in
the record, we print it out:

    import re
    for r in readings:
      if re.search('06', r):
        print r
    Baker 1 2010-06-24      1122.7

So far, this does the same thing as `'06' in r`. But if we want to match
`'06'` or `'07'`, regular expressions let us combine the two comparisons
in a single expression:

    import re
    for r in readings:
      if re.search('06|07', r):
        print r
    Baker 1 2010-06-24      1122.7
    Baker 2 2009-07-24      2819.0

The first argument to `re.search` is the pattern we are searching for,
written as a string. The second argument is the data we are searching
in. It's easy to reverse these accidentally, i.e., to put the data first
and the pattern second. This can be hard to track down, so please be
careful.

The vertical bar in the pattern means "or". It tells regular expression
library that we want to match either the text on the left, or the text
on the right. As we will see [later](#s:mechanics), the regular
expression library can look for both patterns in a single operation.

We are going to be throwing a lot of regular expressions against our
data, so let's write a function that will tell us which records match a
particular pattern. Our function `show_matches` takes a pattern and a
list of strings as arguments. It prints out two stars as a marker if the
pattern matches a string, and just indents with blanks if it does not:

    def show_matches(pattern, strings):
      for s in strings:
        if re.search(pattern, s):
          print '**', s
        else:
          print '  ', s

If we use this function to match `'06|07'` against the data we read in
earlier, it prints stars beside the two records that have month `'06'`
or month `'07'`:

    show_matches('06|07', readings)
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

But if we change the pattern `'06|7'` (without a '0' in front of the
'7'), the pattern seems to match a lot of things that don't have the
month `'06'` or `'07'`:

    show_matches('06|7', readings)
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

To understand why, think back to mathematics. The expression *ab+c*
means "a times b plus c" because multiplication has higher precedence
than addition. If we want to force the other meaning, we have to use
parentheses and write *a(b+c)*.

The same is true for regular expressions. Adjacency has higher
precedence than "or", so the pattern `'06|7'` means, "Either `'06'` or
the digit `'7'`". If we look back at our data, there are a lot of 7's in
our file, and this pattern is matching all of them.

If we want to match `'06'` or `'07'` without repeating the digit '0', we
have to parenthesize it as `'0(6|7)'`. Having said that, most people
probably find the expression `'06|07'` more readable anyway.

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

Matching is useful, but what we really want to do is extract the year,
the month, and the day from our data so that we can reformat them.
Parentheses can help here too: when a regular expression matches a piece
of text, the library automatically remembers what matched against every
parenthesized sub-expression.

Here's a simple example:

    match = re.search('(2009|2010|2011)',
                      'Baker 1\t2009-11-17\t1223.0')
    print match.group(1)

The first string is our pattern. It will match 2009, 2010, or 2011, and
the parentheses around it will make the library remember which of those
three strings was matched. The second string is just the first record
from our data. (Remember, `'\t'` represents a tab.)

When `re.search` is called, it returns `None` if it doesn't find a
match, or a special [match object](glossary.html#match-object) if it
did. The expression `match.group` returns the text that matched the
sub-expression inside the specified set of parentheses counting from the
left. Since this pattern only has one set of parentheses,
`match.group(1)` returns whatever matched what's inside them.

The way sub-expressions are numbered sometimes trips people up. While
Python normally counts from 0, the first match in a regular expression
is extracted with `match.group(1)`, the second with 2, and so forth. The
reason is that `match.group(0)` returns all of the text that the entire
pattern matched.

What if we want to match the month as well as the year? A regular
expression to match legal months would be
`'(01|02|03|04|05|06|07|08|09|10|11|12)'`. An expression to match days
would be three times longer, which would be hard to type and (more
importantly) hard to read.

Instead, we can use the dot character `'.'` to match any single
character. For example, the expression `'....-..-..'` matches exactly
four characters, and `'....-..-..'` matches four characters, a dash, two
more characters, another dash, and two more characters. If we put each
set of dots in parentheses as `'(....)-(..)-(..)'` the three groups
should record the year, month, and day each time there's a successful
match.

Let's test that out by calling `re.search` with the pattern we just
described and the first record from our data:

    match = re.search('(....)-(..)-(..)',
                      'Baker 1\t2009-11-17\t1223.0')
    print match.group(1), match.group(2), match.group(3)
    2009 11 17

When we print out the three groups, we get `'2009'`, `'11'`, and `'17'`,
just as we wanted. Try doing *that* with substring searches…

To recapitulate, leters and digits in a pattern match against
themselves, so `'A'` matches an upper-case A. The vertical bar `'|'`
means "or", a dot `'.'` matches any single character, and we use
parentheses to enforce grouping and to remember things.

Stepping back from the syntax, we have also seen that the right way to
build a pattern is to start with something simple that matches part of
the data we're working with, then add to it piece by piece. We test it
against our data each time we make a change, but also test that it
*doesn't* match things that it shouldn't, because false positive can be
very hard to track down.
