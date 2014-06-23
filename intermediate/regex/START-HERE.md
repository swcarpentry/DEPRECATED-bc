1.  [Simple Patterns](#s:basic)
2.  [Operators](#s:operators)
3.  [Under the Hood](#s:mechanics)
4.  [More Patterns](#s:patterns)
5.  [More Tools](#s:tools)
6.  [One Last Wrinkle](#s:declarative)
7.  [Summing Up](#s:summary)

A couple of years after the Death Star exploded, a hot-shot reporter at
the Daily Planet heard that children in the Shire were starting to act
strangely. Our supervisor sent some grad students off to find out what
was going on. Things didn't go so well for them, but their notebooks
were recovered and later transcribed.

Our job is to read 20 or 30 files, each of which contains several
hundred measurements of background evil levels, and convert them into a
uniform format for further processing. Each of the readings has the name
of the site where the reading was taken, the date the reading was taken
on, and of course the background evil level in millivaders. The problem
is, these files are formatted in different ways. Here is the first one:

~~~~ {fixme="replace with diagram showing characters"}
Site    Date    Evil (millivaders)
----    ----    ------------------
Baker 1 2009-11-17      1223.0
Baker 1 2010-06-24      1122.7
Baker 2 2009-07-24      2819.0
Baker 2 2010-08-25      2971.6
Baker 1 2011-01-05      1410.0
Baker 2 2010-09-04      4671.6
⋮        ⋮               ⋮
~~~~

A single tab character divides the fields in each row into columns. The
site names contain spaces, and the dates are in international standard
format: four digits for the year, two for the month, and two for the
day.

Tabs vs. Spaces

Explain tabs and spaces.

Let's have a look at the second notebook:

~~~~ {fixme="replace with diagram showing characters"}
Site/Date/Evil
Davison/May 22, 2010/1721.3
Davison/May 23, 2010/1724.7
Pertwee/May 24, 2010/2103.8
Davison/June 19, 2010/1731.9
Davison/July 6, 2010/2010.7
Pertwee/Aug 4, 2010/1731.3
Pertwee/Sept 3, 2010/4981.0
⋮        ⋮            ⋮
~~~~

It uses slashes as separators. There don't appear to be spaces in the
site names, but the month names and day numbers vary in length. What's
worse, the months are text, and the order is month-day-year rather than
year-month-day.

We could parse these files using basic string operations, but it would
be difficult. A better approach is to use [regular
expressions](glossary.html#regular-expression). A regular expression is
just a pattern that can match a string. They are actually very common:
when we say `*.txt` to a computer, we mean, "Match all of the filenames
that end in `.txt`." The `*` is a regular expression: it matches any
number of characters.

The rest of this chapter will look at what regular expressions can do,
and how we can use them to handle our data. A warning before we go any
further, though: the notation for regular expressions is ugly, even by
the standards of programming. Se're writing patterns to match strings,
but we're writing those patterns *as* strings using only the symbols
that are on the keyboard, instead of inventing new symbols the way
mathematicians do. The good news is that regular expressions work more
or less the same way in almost every programming language. We will
present examples in Python, but the ideas and notation transfer directly
to Perl, Java, MATLAB, C\#, and Fortran.

Simple Patterns
---------------

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

Operators
---------

Let's go back to those measurements. Notebook \#1 has the site, date,
and background evil level with single tabs as separators. Some of the
site names have spaces, and the dates are in the international standard
format YYYY-MM-DD. However, the fields in Notebook \#2 are separated by
slashes, and use months' names instead of numbers. What's more, some of
the month names are three characters long, while others are four, and
the days are either one or two digits.

Before looking at how to use regular expressions to extract data from
Notebook \#2, let's see how we would do it with simple string
operations. If our records look like `'Davison/May 22, 2010/1721.3'`, we
can split on slashes to separate the site, date, and reading. We could
then split the middle field on spaces to get the month, day, and year,
and then remove the comma from the day if it is present (because some of
our readings don't have a comma after the day).

This is a [procedural](glossary.html#procedural-programming) way to
solve the problem: we tell the computer what procedure to follow step by
step to get an answer. In contrast, regular expressions are
[declarative](glossary.html#declarative-programming): we declare, "This
is what we want," and let the computer figure out how to calculate it.

Our first attempt to parse this data will rely on the `*` operator. It
is a [postfix](glossary.html#postfix-operator) operator, just like the 2
in x^2^, and means, "Zero or more repetitions of the pattern that comes
before it". For example, `'a*'` matches zero or more 'a' characters,
while `'.*'` matches any sequence of characters (including the empty
string) because `'.'` matches anything and `'*'` repeats. Note that the
characters matched by `'.*'` do *not* all have to be the same: the rule
is not, "Match a character against the dot, then repeat that match zero
or more times," but rather, "Zero or more times, match any character."

Here's a test of a simple pattern using `'.*'`:

    match = re.search('(.*)/(.*)/(.*)',
                      'Davison/May 22, 2010/1721.3')
    print match.group(1)
    print match.group(2)
    print match.group(3)

In order for the entire pattern to match, the slashes '/' have to line
up exactly, because '/' only matches against itself. That constraint
ought to make the three uses of `'.*'` match the site name, date, and
reading. Sure enough, the output is:

    Davison
    May 22, 2010
    1271.3

Unfortunately, we've been over-generous. Let's put brackets around each
group in our output to make matches easier to see, then apply this
pattern to the string `'//'`:

    match = re.search('(.*)/(.*)/(.*)',
                      '//')
    print '[' + match.group(1) + ']'
    print '[' + match.group(2) + ']'
    print '[' + match.group(3) + ']'
    []
    []
    []

We don't want our pattern to match invalid records like this (remember,
"Fail early, fail often"). However, `'.*'` can match the empty string
because it is zero occurrences of a character.

Let's try a variation that uses `+` instead of `*`. `+` is also a
postfix operator, but it means "one or more", i.e., it has to match at
least one occurrence of the pattern that comes before it.

    match = re.search('(.+)/(.+)/(.+)',
                      '//')
    print match
    None

As we can see, the pattern `(.+)/(.+)/(.+)` *doesn't* match a string
containing only slashes because there aren't characters before, between,
or after the slashes. And if we go back and check it against valid data,
it seems to do the right thing:

    print re.search('(.+)/(.+)/(.+)',
                    'Davison/May 22, 2010/1721.3')
    print '[' + m.group(1) + ']'
    print '[' + m.group(2) + ']'
    print '[' + m.group(3) + ']'
    [Davison]
    [May 22, 2010]
    [1721.3]

We're going to match a lot of patterns against a lot of strings, so
let's write a function to apply a pattern to a piece of text, report
whether it matches or not, and print out the match groups if it does:

    def show_groups(pattern, text):
      m = re.search(pattern, text)
      if m is None:
        print 'NO MATCH'
        return
      for i in range(1, 1 + len(m.groups())):
        print '%2d: %s' % (i, m.group(i))

We'll test our function against the two records we were just using:

    show_groups('(.+)/(.+)/(.+)',
                'Davison/May 22, 2010/1721.3')
    1: Davison
    2: May 22, 2010
    3: 1721.3

    show_groups('(.+)/(.+)/(.+)',
                '//)
    NO MATCH

All right: if we're using regular expressions to extract the site, date,
and reading, why not add more groups to break up the date while we're at
it?

    show_groups('(.+)/(.+) (.+), (.+)/(.+)',
                'Davison/May 22, 2010/1721.3')
    1: Davison
    2: May
    3: 22
    4: 2010
    5: 1721.3

But wait a second: why doesn't this work?

    show_groups('(.+)/(.+) (.+), (.+)/(.+)',
                'Davison/May 22 2010/1721.3')
    None

The problem is that the string we're trying to match doesn't have a
comma after the day. There is one in the pattern, so matching fails.

We could try to fix this by putting `'*'` after the comma in the
pattern, but that would match any number of consecutive commas in the
data, which we don't want either. Instead, let's use a question mark
`'?'`, which is yet another postfix operator meaning, "0 or 1 of
whatever comes before it". Another way of saying this is that the
pattern that comes before the question mark is optional. If we try our
tests again, we get the right answer in both cases:

    # with comma in data
    show_groups('(.+)/(.+) (.+),? (.+)/(.+)',
                'Davison/May 22, 2010/1721.3')
    1: Davison
    2: May
    3: 22
    4: 2010
    5: 1721.3

    # without comma in data
    show_groups('(.+)/(.+) (.+),? (.+)/(.+)',
                'Davison/May 22 2010/1721.3')
    1: Davison
    2: May
    3: 22
    4: 2010
    5: 1721.3

Let's tighten up our pattern a little bit more. We *don't* want to match
this record:

    Davison/May 22, 201/1721.3

because somebody mis-typed the year, entering three digits instead of
four. (Either that, or whoever took this reading was also using the
physics department's time machine.) We could use four dots in a row to
force the pattern to match exactly four digits:

    (.+)/(.+) (.+),? (....)/(.+)

but this won't win any awards for readability. Instead, let's put the
digit `4` in curly braces `{}` after the dot:

    (.+)/(.+) (.+),? (.{4})/(.+)

In a regular expression, curly braces with a number between them means,
"Match the pattern exactly this many times". Since `.` matches any
character, `.{4}` means "match any four characters".

Let's do a few more tests. Here are some records in which the dates are
either correct or mangled:

    tests = (
        'Davison/May , 2010/1721.3',
        'Davison/May 2, 2010/1721.3',
        'Davison/May 22, 2010/1721.3',
        'Davison/May 222, 2010/1721.3',
        'Davison/May 2, 201/1721.3',
        'Davison/ 22, 2010/1721.3',
        '/May 22, 2010/1721.3',
        'Davison/May 22, 2010/'
    )

And here's a pattern that should match all the records that are correct,
but should fail to match all the records that have been mangled:

    pattern = '(.+)/(.+) (.{1,2}),? (.{4})/(.+)'

We are expecting four digits for the year, and we are allowing 1 or 2
digits for the day, since the expression `{M,N}` matches a pattern from
M to N times.

When we run this pattern against our test data, three records match:

    show_matches(pattern, tests)
    ** Davison/May , 2010/1721.3
    ** Davison/May 2, 2010/1721.3
    ** Davison/May 22, 2010/1721.3
       Davison/May 222, 2010/1721.3
       Davison/May 2, 201/1721.3
       Davison/ 22, 2010/1721.3
       /May 22, 2010/1721.3
       Davison/May 22, 2010/

The second and third matches make sense: 'May 2' and 'May 22' are both
valid. But why does 'May' with no date at all match this pattern? Let's
look at that test case more closely:

    show_groups('(.+)/(.+) (.{1,2}),? (.{4})/(.+)',
                'Davison/May , 2010/1721.3')
    1: Davison
    2: May
    3: ,
    4: 2010
    5: 1721.3

The groups are 'Davison' (that looks right), 'May' (ditto), a ',' on its
own (which is clearly wrong), and then the right year and the right
reading.

Here's what's happened. The space ' ' after 'May' matches the space ' '
in the pattern. The expression "1 or 2 occurrences of any character"
matches the comma ',' because ',' is a character and it occurs once. The
expression ', ' is then not matched against anything, because it's
allowed to match zero characters. '?' means "optional", and in this
case, the regular expression pattern matcher is deciding not to match it
against anything, because that's the only way to get the whole pattern
to match the whole string. After that, the second space matches the
second space in our data. This is obviously not what we want, so let's
modify our pattern again:

    show_groups('(.+)/(.+) ([0-9]{1,2}),? (.{4})/(.+)',
                'Davison/May , 2010/1721.3')
    None

    show_groups('(.+)/(.+) ([0-9]{1,2}),? (.{4})/(.+)',
                'Davison/May 22, 2010/1721.3')
    1: Davison
    2: May
    3: 22
    4: 2010
    5: 1721.3

The pattern `'(.+)/(.+) ([0-9]{1,2}),? (.{4})/(.+)'` does the right
thing for the case where there is no day, and also for the case where
there is one. It works because we have used `[0-9]` instead of `'.'`.

In regular expressions, square brackets `[]` are used to create sets of
characters. For example, the expression `[aeiou]` matches exactly one
vowel, i.e., exactly one occurrence of any character in the set. We can
either write these sets out character by character, as we've done with
vowels, or as "first character '-' last character" if the characters are
in a contiguous range. This is why `'[0-9]'` matches exactly one digit.

Here's our completed pattern:

    (.+)/([A-Z][a-z]+) ([0-9]{1,2}),? ([0-9]{4})/(.+)'

We have added one more feature to it: the name of the month has to begin
with an upper-case letter, i.e., a character in the set `[A-Z]`, which
must followed by one or more lower-case characters in the set `[a-z]`.

This pattern still isn't perfect: the day is one or more occurrences of
the digits 0 through 9, which will allow "days" like '0', '00', and
'99'. It's easiest to check for mistakes like this after we convert the
day to an integer, since trying to handle things like leap years with
regular expressions would be like trying to build a house with a Swiss
army knife.

Finally, the year in our final pattern is exactly four digits, so it's
the set of characters `[0-9]` repeated four times. Again, we will check
for invalid values like '0000' after we convert to integer.

Using the tools we've seen so far, we can write a simple function that
will extract the date from either of the notebooks we have seen so far
and return the year, the month, and the day as strings:

    def get_date(record):
      '''Return (Y, M, D) as strings, or None.'''

      # 2010-01-01
      m = re.search('([0-9]{4})-([0-9]{2})-([0-9]{2})',
                    record)
      if m:
        return m.group(1), m.group(2), m.group(3)

      # Jan 1, 2010 (comma optional, day may be 1 or 2 digits)
      m = re.search('/([A-Z][a-z]+) ([0-9]{1,2}),? ([0-9]{4})/',
                    record)
      if m:
        return m.group(3), m.group(1), m.group(2)

      return None

We start by testing whether the record contains an ISO-formatted date
YYYY-MM-DD. If it does, then we return those three fields right away.
Otherwise, we test the record against a second pattern to see if we can
find the name of a month, one or two digits for the day, and four digits
for the year with slashes between the fields. If so, we return what we
find, permuting the order to year, month, day. Finally, if neither
pattern matched we return `None` to signal that we couldn't find
anything in the data.

This is probably the most common way to use regular expressions: rather
than trying to combine everything into one enormous pattern, we have one
pattern for each valid case. We test of those cases in turn; if it
matches, we return what we found, and if it doesn't, we move on to the
next pattern. Writing our code this way make it easier to understand
than using a single monster pattern, and easier to extend if we have to
handle more data formats.

Under the Hood
--------------

The regular expression `'([A-Z][a-z]+) ([0-9]{1,2}),? ([0-9]{4})'`
matches a single upper-case character and one or more lower-case
characters, a space, one or two digits, an optional comma, another
space, and exactly four digits. That is pretty complex, and knowing a
little about how the computer actually does it will help us debug
regular expressions when they don't do what we want.

Regular expressions are implemented using [finite state
machines](glossary.html#finite-state-machine). Here's a very simple FSM
that matches exactly one lower case 'a':

![FSM matching a single lower case
'a'](../img/regexp/fsm-single-lower-case-a.png)

Matching starts with the incoming arrow on the left, which takes us to
the first state in our finite state machine. The only way to get from
there to the second state is to match the 'a' on the arc between states
1 and 2. The dot in the middle of the second state means that it's an
end state. We must be in one of these states at the end of our match in
order for the match to be valid.

Now that we have an FSM that matches the very simple regular expression
`'a'`, let's see if we can do something a little more interesting.
Here's a finite state machine that matches one or more occurrences of
the letter 'a':

![FSM matching one or more letter
'a'](../img/regexp/fsm-one-or-more-a.png)

The first arc labelled 'a' gets us from the initial state to an end
state, but we don't have to stop there: the curved arc at the top allows
us to match another 'a', and brings us back to the same state. We can
then match another 'a', and another, and so on indefinitely. (Note that
we don't have to stop in the end state the first time we reach it: we
just have to be in an end state when we run out of input.) The pattern
this FSM matches is `'a+'`, since one 'a' followed by zero or more
others is the same as one or more occurences of 'a'.

Here's another FSM that matches against the letter 'a' or nothing:

![FSM matching one letter 'a' or
nothing](../img/regexp/fsm-one-a-or-nothing.png)

The top arc isn't marked, so that transition is free: we can go from the
first state to the second state without consuming any of our input. This
is "a or nothing", which is the same as `'a?'`, i.e., an optional
character 'a'.

This regular expression looks like the one that matches 'a' one or more
times, except there is an extra arc to get us from the first state to
the second without consuming any input:

![FSM matching zero or more letter
'a'](../img/regexp/fsm-zero-or-more-a.png)

It is therefore equivalent to the pattern `'a*'`, i.e., it matches
nothing at all (taking that free transition from the first state to the
second) or one or more occurrences of 'a'. We can simplify this
considerably like this:

![FSM matching zero or more letter
'a'](../img/regexp/fsm-simpler-zero-or-more-a.png)

The simple FSMs we have seen so far are enough to implement most of the
regular expressions in the previous sections. For example, look at this
finite state machine:

![A more complex FSM](../img/regexp/fsm-complex.png)

We can either take the top route or the bottom. The top route is `a+`;
the bottom route is a 'b', followed by either a 'c' or a 'd', so this
whole thing is equivalent to the regular expression `'a+|(b(c|d))'`. An
input string that matches any of these paths will leave us in that final
end state.

The most important thing about finite state machines is that the action
they take at a node depends on only the arcs out of that node and the
characters in the target data. Finite state machines do *not* remember
how they got to a particular node: decision-making is always purely
local.

This means that there are many patterns that regular expressions
*cannot* match. For example, it is impossible to write a regular
expression to check if nested parentheses match. If we want to see
whether '(((…)))' is balanced, we need some kind of memory, or at least
a counter, and there isn't any in a finite state machine.

Similarly, if we want to check whether a word contains each vowel only
once, the only way to do it is to write a regular expression with 120
clauses, that checks for each possible permutation of 'aeiou'
explicitly. We cannot write a regular expression that matches an
arbitrary vowel, and then subtracts that vowel from the set of vowels
yet to be matched, because once again, that would require some kind of
memory, and finite state machines don't have any.

Despite this limitation, regular expressions are tremendously useful.
The first reason is that they are really fast. After the computer does
some pre-calculation (essentially, once it turns the regular expression
into a finite state machine) a regular expression can be matched against
input by looking at each input character only once. That means that the
time required to find patterns with regular expressions grows in
proportion to the size of the data. The time required for most other
pattern-matching techniques grows much faster, so if regular expressions
can do the job, they are almost always the most efficient option
available.

Another reason for using regular expressions is that they are more
readable than other alternatives. You might not think so looking at the
examples so far, but imagine writing lines of code to match that same
patterns. Nobody would claim that regular expressions are easy to
understand, but they're a lot easier than two dozen lines of substring
operations.

More Patterns
-------------

Now that we know how regular expressions work, let's have a look at
Notebook \#3:

    Date Site Evil(mvad)
    May 29 2010 (Hartnell) 1029.3
    May 30 2010 (Hartnell) 1119.2
    June 1 2010 (Hartnell) 1319.4
    May 29 2010 (Troughton) 1419.3
    May 30 2010 (Troughton) 1420.0
    June 1 2010 (Troughton) 1419.8
    ⋮            ⋮           ⋮

It has the date as three fields, the site name in parentheses, and then
the reading. We know how to parse dates in this format, and the fields
are separated by spaces, but how do we match those parentheses? The
parentheses we have seen in regular expressions so far haven't matched
characters: they have created groups.

The way we solve this problem—i.e., the way we match a literal left
parenthesis '(' or right parenthesis ')'—is to put a backslash in front
of it. This is another example of an [escape
sequence](glossary.html#escape-sequence): just as we use the
two-character sequence `'\t'` in a string to represent a literal tab
character, we use the two-character sequence `'\('` or `'\)'` in a
regular expression to match the literal character '(' or ')'.

To get that backslash '\\' into the string, though, we have to escape
*it* by doubling it up. This has nothing to do with regular expressions:
it is Python's rule for putting backslashes in strings. The string
representation of the regular expression that matches an opening
parenthesis is therefore `'\\('`. This can be confusing, so let's take a
look at the various layers involved.

Our program text—i.e., what's stored in our `.py` file—looks like this:

    # find '()' in text
    m = re.search('\\(\\)', text)
    ⋮    ⋮    ⋮

The file has two backslashes, an open parenthesis, two backslashes, and
a close parenthesis inside quotes:

  ---- ---- --- ---- ---- ---
  \\   \\   (   \\   \\   )
  ---- ---- --- ---- ---- ---

When Python reads that file in, it turns the two-character sequence
'\\\\' into a single '\\' character in the string in memory. This is the
first round of escaping.

  ---- --- ---- ---
  \\   (   \\   )
  ---- --- ---- ---

When we hand that string '\\(\\)' to the regular expression library, it
takes the two-character sequence '\\(' and turns it into an arc in the
finite state machine that matches a literal parenthesis:

![More complex FSM](../img/regexp/fsm-match-parentheses.png)

Turning this over, if we want a literal parenthesis to be matched, we
have to give the regular expression library '\\('. If we want to put
'\\(' in a string, we have to write it in our `.py` file as '\\\\('.

With that out of the way, let's go back to Notebook \#3. The regular
expression that will extract the five fields from each record is
`'([A-Z][a-z]+) ([0-9]{1,2}) ([0-9]{4}) \\((.+)\\) (.+)'`, which is:

-   a word beginning with an upper-case character followed by one or
    more lower-case characters,
-   a space,
-   one or two digits,
-   another space,
-   four digits,
-   another space,
-   some stuff involving backslashes and parentheses,
-   another space,
-   and then one or more characters making up the reading.

If we take a closer look at that "stuff", `'\\('` and `'\\)'` are how we
write the regular expressions that match a literal open parenthesis '('
or close parenthesis ')' character in our data. The two inner
parentheses that don't have backslashes in front of them create a group,
but don't match any literal characters. We create that group so that we
can save the results of the match (in this case, the name of the site).

Now that we know how to work with backslahes in regular expressions, we
can take a look at character sets that come up frequently enough to
deserve their own abbreviations. If we use `'\d'` in a regular
expression it matches the digits 0 through 9. If we use `'\s'`, it
matches the whitespace characters, space, tab, carriage return, and
newline. `'\w'` matches word characters; it's equivalent to the set
`'[A-Za-z0-9_]'` of upper-case letters, lower-case letters, digits, and
the underscore. (It's actually the set of characters that can appear in
a variable name in a programming language like C or Python.) Again, in
order to write one of these regular expressions as a string in Python,
we have to double the backslashes.

Now that we've seen these character sets, we can take a look at an
example of really bad design. The regular expression `'\S'` means
"non-space characters", i.e., everything that *isn't* a space, tab,
carriage return, or newline. That might seem to contradict what we said
in the previous paragraph, but if we look closely, that's an upper-case
'S', not a lower-case 's'.

Similarly, and equally unfortunately, `'\W'` means "non-word characters"
provided it's an upper-case 'W'. Upper- and lower-case 'S' and 'W' look
very similar, particularly when there aren't other characters right next
to them to give context. This means that these sequences are very easy
to mis-type and mis-read. Everyone eventually uses an upper-case 'S'
when they meant to use a lower-case 's' or vice versa, and then wastes a
few hours trying to track it down. So please, if you're ever designing a
library that's likely to be widely used, try to choose a notation that
doesn't make mistakes this easy.

Along with the abbreviations for character sets, the regular expression
library recognizes a few shortcuts for things that aren't actual
characters. For example, if we put a circumflex `'^'` at the start of a
pattern, it matches the beginning of the input text. (Note that there's
no backslash in front of it.) This means that the pattern `'^mask'` will
match the text `'mask size'`, because the letters 'mask' come at the
start of the string, but will *not* match the word `'unmask'`. Going to
the other end, if dollar sign `'$'` is the last character in the
pattern, it matches the end of the line rather than a literal '\$', so
'temp\$' will match the string 'high-temp', but not the string
'temperature'.

Regular Expressions and Newlines

The full rule is slightly more complicated. By default, regular
expressions act as if newline characters were the ends of records. For
example, the `'.'` pattern matches everything *except* a newline. This
normally doesn't matter, since most I/O routines return one line of text
at a time, but if we read a whole file into a single string, then try to
match across line boundaries, we may not get the behavior we expect. We
can use the `MULTILINE` option in our matches to prevent this; please
see the regular expression documentation for details.

A third shortcut that's often useful is `'\b'`, often called "break". It
doesn't match any characters; instead, it matches the boundary between
word and non-word characters (where "word" means upper and lower case
characters, digits, and the underscore). For example, the regular
expression `'\bage\b'` will match the string `'the age of'` because
there's a non-word character right before the 'a' and another non-word
character right after the 'e'. That same pattern will not match the word
`'phage'` because there isn't a transition from non-word to word
characters, or vice versa, right before the 'a'. And remember: to get
that regular expression int our program, we have to escape the
backslashes using `'\\bage\\b'`.

One Last Wrinkle
----------------

Let's have one last look at the function we wrote to extract data from
lab notebooks:

    def get_date(record):
      '''Return (Y, M, D) as strings, or None.'''

      # 2010-01-01
      m = re.search('([0-9]{4})-([0-9]{2})-([0-9]{2})',
                    record)
      if m:
        return m.group(1), m.group(2), m.group(3)

      # Jan 1, 2010 (comma optional, day may be 1 or 2 digits)
      m = re.search('/([A-Z][a-z]+) ([0-9]{1,2}),? ([0-9]{4})/',
                    record)
      if m:
        return m.group(3), m.group(1), m.group(2)

      return None

We can make it easier to add new patterns to this function by making it
more declarative. The trick is to combine the regular expressions and
the IDs of the groups we want to return:

    def get_fields(record):
      '''Return (Y, M, D, site, reading) or None.'''

      patterns = [
        ['(.+)\t([0-9]{4})-([0-9]{2})-([0-9]{2})\t(.+)',      2, 3, 4, 1, 5],
        ['(.+)/([A-Z][a-z]+) ([0-9]{1,2}),? ([0-9]{4})/(.+)', 4, 2, 3, 1, 5]
      ]
      for pattern, year, month, day, site, reading in patterns:
        m = re.search(pattern, record)
        if m:
          return m.group(year), m.group(month), m.group(day),
                 m.group(site), m.group(reading)

      return None

Each entry in the list `patterns` has two parts: a regular expression,
and then the indices of the group that will contain the year, month,
day, site, and reading if that pattern matches. The loop tries the
regular expressions in `patterns` one by one. As soon as a pattern
matches it returns the matched groups, permuting them according to the
indices so that the data always comes back in the same order. To handle
the format in Notebook \#3, we just add one line to this table:

        ['([A-Z][a-z]+) ([0-9]{1,2}) ([0-9]{4}) \\((.+)\\) (.+)', 3, 1, 2, 4, 5]

Using a table might not seem like much of an improvement over the
"match, extract, and return" style we have been using so far. However,
the table-based approach has one major advantage: it signals to the
reader that all the patterns are being handled the same way. It's all
too easy for programmers to tweak the branches of a "match, extract, and
return" function so that each possibility is handled in a slightly
different way. This makes it very hard for readers to understand what's
going on, and equally hard for the next programmer in line to debug or
extend the code. The more declarative the code is, the more confidence
readers can have that there really is only one thing for them to
understand.

More Tools
----------

To end our exploration of regular expressions, let's work through a
moderately complex problem and introduce a few more tools in the regular
expression library. Our starting point is an archive of several thousand
papers and theses written in LaTeX, a text-based document formatting
program. LaTeX documents use labels to refer to items in a shared
bibliography. Our job is to find out how often citations appear
together, i.e., how often paper X is cited in the same document as paper
Y. To answer this question we need to extract the set of citation labels
from each document.

Let's have a closer look at our input:

    Granger's work on graphs \cite{dd-gr2007,gr2009},
    particularly ones obeying Snape's Inequality
    \cite{ snape87 } (but see \cite{quirrell89}),
    has opened up new lines of research.  However,
    studies at Unseen University \cite{stibbons2002,
    stibbons2008} highlight several dangers.
    ⋮    ⋮    ⋮

Citations in LaTeX are written using `\cite{…}`, with cross-reference
labels in the curly braces. A single citation can include two or more
labels separated by commas. There may be white space before or after
labels or line breaks where a citation is split across two lines, and
there can be multiple citations per line.

Our first idea is to use a group to capture everything inside the curly
braces following the word 'cite':

    m = re.search('cite{(.+)}', 'a \\cite{X} b')
    print m.groups()

~~~~ {.out}
('X',)
~~~~

It seems to work in one simple case, but what if there are multiple
citations on a single line?

    m = re.search('cite{(.+)}', 'a \\cite{X} b \\cite{Y} c')
    print m.groups()
    ('X} b \\cite{Y',)

It looks like we're capturing the text *between* the citations. The
reason is that regular expression matching is
[greedy](glossary.html#greedy-matching): it matches as much text as it
can, and the '.' in `'.+'` will match all the characters from the first
curly brace to the last one, including the intervening citations and
curly braces.

The diagnosis of our problem suggests its solution: let's have the
regular expression match everything *except* a closing curly brace. This
is easy to do: if the first character of a set in square brackets is the
circumflex '\^', then the set is negated, i.e., it matches everything
*except* the characters in the set. The expression `[^}]` therefore
matches every character except a closing curly brace. Let's try it out:

    m = re.search('cite{([^}]+)}', 'a \\cite{X} b')
    print m.groups()
    ('X,)

This works for a single citation: all we've done is change '.' to the
negated set. What about multiple citations on a single line?

    m = re.search('cite{([^}]+)}', 'a \\cite{X} b \\cite{Y} c')
    print m.groups()
    ('X,)

It's not gobbling up text we don't want it to, but it's only capturing
the first citation. Somehow, we need to extract all matches, not just
the first.

The regular expression library has a function to do exactly this: if we
use `re.findall` instead of `re.search`, it will give us back a list of
all the substrings that matched our pattern. Remember, whatever your
problem is, someone has probably run into it before, and there's
probably something in the library to help you. Knowing what's in the
library is as important to a programmer as knowing what's in the
literature is to a scientist. The bad news is, it's usually hard to find
things in libraries or their documentation unless you already know
enough about your problem to know what keywords to search for.

Let's give `findall` a try:

    print re.findall('cite{([^}]+)}', 'a \\cite{X} b \\cite{Y} c')
    ['X', 'Y']

It seems to produce the right output—not bad for a 7-character change.
What about spaces in citations?

    print re.search('cite{([^}]+)}', 'a \\cite{ X} b \\cite{Y } c').groups()
    [' X', 'Y ']

The good news is, nothing breaks. The bad news is, the spaces are saved
by `findall`, which isn't really what we want. We could tidy this up
after the fact using `string.strip`, but let's modify the pattern
instead:

    print re.findall('cite{\\s*([^}]+)\\s*}', 'a \\cite{ X} b \\cite{Y } c')
    ['X', 'Y ']

If you recall, `'\s'` is an abbrevation for the set of whitespace
characters, so these uses of `'\s*'` match zero or more spaces
immediately after the opening curly brace, or immediately before the
closing one (and as always we have to write `'\\s'` to get the backslash
into the Python string). However, the space after the 'Y' is still being
returned in the matched text.

Once again, the problem is that regular expressions are greedy: the
space after the 'Y' isn't a closing curly brace, so it's matched by the
negated character set, and included in the returned string. The `'\s*'`
that's supposed to match the trailing space is then matched against zero
characters instead of one. It's not what we want, but it's legal.

Let's force our match to line up with the break from word to non-word
characters using `'\b'`:

    print re.findall('cite{\\s*\\b([^}]+)\\b\\s*}', 'a \\cite{ X} b \\cite{Y } c')
    ['X', 'Y']

It works! check this last example: in the PowerPoint, there's still a
space before the 'X' The change is to put `'\b'` after the first
unwanted spaces, and before the last ones. Since the curly braces around
the citation labels are also non-word characters, the pattern matches
even if there aren't any opening or trailing spaces.

The last hurdle is to handle multiple labels inside a single pair of
curly braces. The pattern we've built so far doesn't explode when there
are two or more labels, and it even handles spaces after the commas, but
it returns all those labels as a single lump of text:

    print re.findall('cite{\\s*\\b([^}]+)\\b\\s*}', '\\cite{X,Y} ')
    ['X,Y']

    print re.findall('cite{\\s*\\b([^}]+)\\b\\s*}', '\\cite{X, Y, Z} ')
    ['X, Y, Z']

We actually could write a pattern that would break everything up on
commas, but it would need some very advanced features of the regular
expression library. Instead, let's use another basic function,
`re.split`, to separate multiple labels. `re.split` does the same thing
as `string.split`, but unlike its simpler cousin, it breaks things up
everywhere that a pattern matches.

The best way to show how it works is to write the function we originally
set out to create. Let's start with a skeleton that includes some test
data, a function that does nothing (but doesn't just fail), and a couple
of lines that call that function and display the result:

    def get_citations(text):
      '''Return the set of all citation tags found in a block of text.'''
      return set() # to be done

    if __name__ == '__main__':
      test = '''\
    Granger's work on graphs \cite{dd-gr2007,gr2009},
    particularly ones obeying Snape's Inequality
    \cite{ snape87 } (but see \cite{quirrell89}),
    has opened up new lines of research.  However,
    studies at Unseen University \cite{stibbons2002,                                                                                                       
    stibbons2008} highlight several dangers.'''

      print get_citations(test)
    set([])

Now let's write our function. For readability's sake, we'll put our
patterns at the top and give them memorable names. Inside the function,
we'll pull out all the citations using the first pattern, then split
each result everywhere there's a comma with optional spaces before or
after it. We'll stuff all the results into a set, and return that. If no
matches were found, that set will be empty.

    import re

    p_cite = 'cite{\\s*\\b([^}]+)\\b\\s*}'
    p_split = '\\s*,\\s*'

    def get_citations(text):
      '''Return the set of all citation tags found in a block of text.'''

      result = set()
      match = re.findall(p_cite, text)
      if match:
        for citation in match:
          cites = re.split(p_split, citation)
          for c in cites:
            result.add(c)

      return result

We can use one more trick from the regular expression library to make
this function more efficient. Instead of turning the regular expression
into a finite state machine over and over again, we can compile the
regular expression and save the resulting object:

    import re

    p_cite = re.compile('cite{\\s*\\b([^}]+)\\b\\s*}')
    p_split = re.compile('\\s*,\\s*')

    def get_citations(text):
      '''Return the set of all citation tags found in a block of text.'''

      result = set()
      match = p_cite.findall(text)
      if match:
        for citations in match:
          label_list = p_split.split(citations)
          for label in label_list:
            result.add(label)

      return result

That object has methods with the same names as the functions we've been
using from the library, like `search` and `findall`, but if we're using
the same pattern over and over again, compiling it once and re-using the
compiled object is much faster.

As you can see, the changes required are very small: instead of saving
the textual representations of our expressions, we compile them, and
then instead of calling the top-level functions from the regular
expression library, we call the methods of those saved objects. The
result is a set of all the citations in our test data, pulled out with
just a dozen lines of code:

    if __name__ == '__main__':
      test = '''\
    Granger's work on graphs \cite{dd-gr2007,gr2009},
    particularly ones obeying Snape's Inequality
    \cite{ snape87 } (but see \cite{quirrell89}),
    has opened up new lines of research.  However,
    studies at Unseen University \cite{stibbons2002,                                                                                                       
    stibbons2008} highlight several dangers.'''

      print get_citations(test)
    set(['gr2009', 'stibbons2002', 'dd-gr2007', 'stibbons2008',
           'snape87', 'quirrell89'])

Finally, if we are going to compile our regular expressions, we can make
them even easier to understand by using *verbose mode* to add comments.
Verbose mode tells Python to ignore whitespace and comments in the
regular expression, which lets us write patterns like this:

    p_cite = '''
        cite{          # start with literal 'cite{'
        \\s*           # then some optional spaces
        \\b            # up to a start-of-word boundary
        ([^}]+)        # then anything that isn't a closing '}'
        \\b            # then an end-of-word boundary
        \\s*           # and some more optional spaces
        }              # and the closing '}'
    '''
    matcher = re.compile(p_cite, re.VERBOSE)

Documenting patterns like this makes them much easier to fix and extend.

Summing Up
----------

Every regular expressions is actually a little program: a piece of data
that tells a computer how to behave. The steps the computer goes through
when executing a Python or Fortran program are more complicated than the
ones it uses to compile and apply a regular expression, but the
principles are exactly the same.

One last point to take away from this chapter is that if we know we are
going to use regular expressions to read in data, we should choose a
format for that data that's easy for regular expressions to match.
Optional commas, tabs that might be repeated, and other things that make
data easy for people to type in actually make it harder for programs to
read that data reliably. This tension between what's easy for the
machine and what's easy for the user never goes away, but if we're
conscious of it, we can find a happy medium.
