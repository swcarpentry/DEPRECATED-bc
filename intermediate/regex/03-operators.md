---
layout: lesson
root: ../..
title: Operators
---


Let's go back to those measurements. Notebook  1 has the site, date,
changes as suggested by @ctjacobs
and background evil level with single tabs as separators. Some of the
site names have spaces, and the dates are in the international standard
format YYYY-MM-DD. However, the fields in Notebook  2 are separated by
slashes, and use month names instead of numbers. What's more, some of
the month names are three characters long, while others are four, and
the days are either one or two digits.

### Using simple string operations gets tedious quickly

Before looking at how to use regular expressions to extract data from
Notebook 2, let's see how we would do it with simple string
operations. If our records look like `'Davison/May 22, 2010/1721.3'`, we
can split on slashes to separate the site, date, and reading. We could
then split the middle field on spaces to get the month, day, and year,
and then remove the comma from the day if it is present (because some of
our readings don't have a comma after the day).

This is a [procedural](../../gloss.html#procedural-programming) way to
solve the problem: we tell the computer what procedure to follow step by
step to get an answer. In contrast, regular expressions are
[declarative](../../gloss.html#declarative-programming): we declare, "This
is what we want," and let the computer figure out how to calculate it.

Good definitions rely on us being able to define characters that stand in place for others.
Having to explicitly define the literal string usually isn't a big help,
what we want is someway of defining the general pattern.
This is where operators in regular expressions come in handy. 

### Operators specify patterns that simplify regular expressions

Operators are the bread and butter of regular expressions.
They are simply characters that specify other patterns of characters
(sometimes of varying length).

You'll have seen operators in use before.
The `*` operator familiar from many a GUI's find box or command-line wildcards is a very common one.
Many combinations of characters such as `\s` are operators too.

Here is a quick quiz.
What happens if I try to use the pattern `txt/files/(*.txt)`
on this text `txt/files/file.txt`?

1.  It matches `file.txt`
2.  It matches the whole string
3.  It won't work at all (the regex won't compile)

Perhaps surprisingly, the answer is `3`.
The regex won't compile because the `.` and `*` are operators that don't mean the same as they do in GUI search boxes
and are a source of gotcha's when building regular expressions.

### Using operators

Our first attempt to parse this data will rely on the `*` operator.
It is a [postfix](../../gloss.html#postfix-operator) operator,
and means "Zero or more repetitions of the pattern that comes
before it". For example, `a*` matches zero or more `a` characters,
while `.*` matches any sequence of characters (including the empty
string) because `.` matches anything and `*` repeats. Note that the
characters matched by `.*` do *not* all have to be the same: the rule
is not, "Match a character against the dot, then repeat that match zero
or more times," but rather, "Zero or more times, match any character."

Here's a test of a simple pattern using `.*`:


~~~
match = re.search('(.*)/(.*)/(.*)',
                  'Davison/May 22, 2010/1721.3')
print match.group(1)
print match.group(2)
print match.group(3)
~~~
{:class="in"}

In order for the entire pattern to match, the slashes `/` have to line
up exactly, because '/' only matches against itself. That constraint
ought to make the three uses of `.*` match the site name, date, and
reading. Sure enough, the output is:


~~~
Davison
May 22, 2010
1271.3
~~~
{:class="out"}

Unfortunately, we've been over-generous. Let's put brackets around each
group in our output to make matches easier to see, then apply this
pattern to the string `'//'`:

~~~
match = re.search('(.*)/(.*)/(.*)',
                  '//')
print '[' + match.group(1) + ']'
print '[' + match.group(2) + ']'
print '[' + match.group(3) + ']'
~~~
{:class="in"}

~~~
[]
[]
[]
~~~
{:class="out"}

We don't want our pattern to match invalid records like this (remember,
"Fail early, fail often"). However, `'.*'` can match the empty string
because it is zero occurrences of a character.

Let's try a variation that uses `+` instead of `*`. `+` is also a
postfix operator, but it means "one or more", i.e., it has to match at
least one occurrence of the pattern that comes before it.

~~~
match = re.search('(.+)/(.+)/(.+)',
                  '//')
print match
~~~
{:class="in"}

~~~
None
~~~
{:class="out"}

As we can see, the pattern `(.+)/(.+)/(.+)` *doesn't* match a string
containing only slashes because there aren't characters before, between,
or after the slashes. And if we go back and check it against valid data,
it seems to do the right thing:


~~~
print re.search('(.+)/(.+)/(.+)',
                'Davison/May 22, 2010/1721.3')
print '[' + m.group(1) + ']'
print '[' + m.group(2) + ']'
print '[' + m.group(3) + ']'
~~~
{:class="in"}

~~~
[Davison]
[May 22, 2010]
[1721.3]
~~~
{:class="out"}

We're going to match a lot of patterns against a lot of strings, so
let's write a function to apply a pattern to a piece of text, report
whether it matches or not, and print out the match groups if it does:

~~~
def show_groups(pattern, text):
    m = re.search(pattern, text)
    if m is None:
        print 'NO MATCH'
        return
    for i in range(1, 1 + len(m.groups())):
        print '%2d: %s' % (i, m.group(i))
~~~
{:class="in"}


We'll test our function against the two records we were just using:

~~~
show_groups('(.+)/(.+)/(.+)',
            'Davison/May 22, 2010/1721.3')
~~~
{:class="in"}

~~~
1: Davison
2: May 22, 2010
3: 1721.3
~~~
{:class="out"}

~~~
show_groups('(.+)/(.+)/(.+)',
            '//')
~~~
{:class="in"}

~~~
NO MATCH
~~~
{:class="out"}

All right: if we're using regular expressions to extract the site, date,
and reading, why not add more groups to break up the date while we're at
it?

~~~
show_groups('(.+)/(.+) (.+), (.+)/(.+)',
            'Davison/May 22, 2010/1721.3')
~~~
{:class="in"}

~~~
1: Davison
2: May
3: 22
4: 2010
5: 1721.3
~~~
{:class="out"}

But wait a second: why doesn't this work?

~~~
show_groups('(.+)/(.+) (.+), (.+)/(.+)',
            'Davison/May 22 2010/1721.3')
~~~
{:class="in"}

~~~
None
~~~
{:class="out"}

The problem is that the string we're trying to match doesn't have a
comma after the day. There is one in the pattern, so matching fails.

We could try to fix this by putting `*` after the comma in the
pattern, but that would match any number of consecutive commas in the
data, which we don't want either. Instead, let's use a question mark
`?`, which is yet another postfix operator meaning, "0 or 1 of
whatever comes before it". Another way of saying this is that the
pattern that comes before the question mark is optional. If we try our
tests again, we get the right answer in both cases:

~~~
# with comma in data
show_groups('(.+)/(.+) (.+),? (.+)/(.+)',
            'Davison/May 22, 2010/1721.3')
~~~
{:class="in"}

~~~
1: Davison
2: May
3: 22
4: 2010
5: 1721.3
~~~
{:class="out"}

~~~
# without comma in data
show_groups('(.+)/(.+) (.+),? (.+)/(.+)',
            'Davison/May 22 2010/1721.3')
~~~
{:class="in"}

~~~
1: Davison
2: May
3: 22
4: 2010
5: 1721.3
~~~
{:class="out"}

Let's tighten up our pattern a little bit more. We *don't* want to match
this record:

~~~
Davison/May 22, 201/1721.3
~~~

because somebody mis-typed the year, entering three digits instead of
four. (Either that, or whoever took this reading was also using the
physics department's time machine.) We could use four dots in a row to
force the pattern to match exactly four digits:

~~~
(.+)/(.+) (.+),? (....)/(.+)
~~~

but this won't win any awards for readability. Instead, let's put the
digit `4` in curly braces `{}` after the dot:

~~~
(.+)/(.+) (.+),? (.{4})/(.+)
~~~

In a regular expression, curly braces with a number between them means,
"Match the pattern exactly this many times". Since `.` matches any
character, `.{4}` means "match any four characters".

Let's do a few more tests. Here are some records in which the dates are
either correct or mangled:


~~~
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
~~~
{:class="in"}

And here's a pattern that should match all the records that are correct,
but should fail to match all the records that have been mangled:

~~~
pattern = '(.+)/(.+) (.{1,2}),? (.{4})/(.+)'
~~~
{:class="in"}

We are expecting four digits for the year, and we are allowing 1 or 2
digits for the day, since the expression `{M,N}` matches a pattern from
M to N times.

When we run this pattern against our test data, three records match:

~~~
show_matches(pattern, tests)
~~~
{:class="in"}

~~~
** Davison/May , 2010/1721.3
** Davison/May 2, 2010/1721.3
** Davison/May 22, 2010/1721.3
   Davison/May 222, 2010/1721.3
   Davison/May 2, 201/1721.3
   Davison/ 22, 2010/1721.3
   /May 22, 2010/1721.3
   Davison/May 22, 2010/
~~~
{:class="out"}

The second and third matches make sense: `May 2` and `May 22` are both
valid. But why does `May` with no date at all match this pattern? Let's
look at that test case more closely:

~~~
show_groups('(.+)/(.+) (.{1,2}),? (.{4})/(.+)',
            'Davison/May , 2010/1721.3')
~~~
{:class="in"}

~~~
1: Davison
2: May
3: ,
4: 2010
5: 1721.3
~~~
{:class="out"}

The groups are `Davison` (that looks right), `May` (ditto), a `,` on its
own (which is clearly wrong), and then the right year and the right
reading.

Here's what's happened.  The space ` ` after `May` matches the space
` ` in the pattern.  The expression "1 or 2 occurrences of any
character" matches the comma `,` because `,` is a character and it
occurs once. The expression `, ` is then not matched against anything,
because it's allowed to match zero characters. `?` means "optional",
and in this case, the regular expression pattern matcher is deciding
not to match it against anything, because that's the only way to get
the whole pattern to match the whole string. After that, the second
space matches the second space in our data. This is obviously not what
we want, so let's modify our pattern again:

~~~
show_groups('(.+)/(.+) ([0-9]{1,2}),? (.{4})/(.+)',
            'Davison/May , 2010/1721.3')
~~~
{:class="in"}

~~~
None
~~~
{:class="out"}


~~~
show_groups('(.+)/(.+) ([0-9]{1,2}),? (.{4})/(.+)',
            'Davison/May 22, 2010/1721.3')
~~~
{:class="in"}

~~~
1: Davison
2: May
3: 22
4: 2010
5: 1721.3
~~~
{:class="out"}

The pattern `(.+)/(.+) ([0-9]{1,2}),? (.{4})/(.+)` does the right
thing for the case where there is no day, and also for the case where
there is one. It works because we have used `[0-9]` instead of `.`.

In regular expressions, square brackets `[]` are used to create sets of
characters (sometimes called character classes). For example, the expression `[aeiou]` matches exactly one
vowel, i.e., exactly one occurrence of any character in the set. We can
either write these sets out character by character, as we've done with
vowels, or as "first character `-` last character" if the characters are
in a contiguous range. This is why `[0-9]` matches exactly one digit.

Here's our completed pattern:

~~~
(.+)/([A-Z][a-z]+) ([0-9]{1,2}),? ([0-9]{4})/(.+)
~~~

We have added one more feature to it: the name of the month has to begin
with an upper-case letter, i.e., a character in the set `[A-Z]`, which
must followed by one or more lower-case characters in the set `[a-z]`.

This pattern still isn't perfect: the day is one or more occurrences of
the digits `0` through `9`, which will allow "days" like `0`, `00`, and
`99`. It's easiest to check for mistakes like this after we convert the
day to an integer, since trying to handle things like leap years with
regular expressions would be like trying to build a house with a Swiss
army knife.

Finally, the year in our final pattern is exactly four digits, so it's
the set of characters `[0-9]` repeated four times. Again, we will check
for invalid values like `0000` after we convert to integer.

Using the tools we've seen so far, we can write a simple function that
will extract the date from either of the notebooks we have seen so far
and return the year, the month, and the day as strings:


~~~
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
~~~
{:class="in"}

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

### Anchors

It is possible to 'anchor' the pattern to a particular part of the
string, so that it can only match in one region, like at the start or
end of the string. the `^` anchor will match the subsequent pattern
only at the start of a string. Likewise the `$` operator will match
the previous pattern only at the end of a line. Let's look at a
contrived example and imagine that we're only interested in data from
one site.

~~~
m = re.search('(^Davison.*)', 'Davison/May 22, 2010/1721.3')
print m.group(1)
~~~
{:class="in"}
~~~
Davison/May 22, 2010/1721.3
~~~
{:class="out"}

Whereas,

~~~
m = re.search('(^Baker.*)', 'Davison/May 22, 2010/1721.3')
print m
~~~
{:class="in"}
~~~
None
~~~
{:class="out"}

Likewise, if we switch the order of columns:

~~~
m = re.search('(^Davison)', '1721.3/May 22, 2010/Davison')
print m
~~~
{:class="in"}

~~~    
None
~~~
{:class="out"}

Since `^Davison` will only match the occurrence at the beginning of the string. Matching the end of the string behaves similarly:

~~~
m = re.search('(.*Davison$)', '1721.3/May 22, 2010/Davison')
print m.group(1)
~~~
{:class="in"}

### Metacharacters

A common thing in regular expressions are metacharacters. These are special pairs of characters that denote classes of characters and help make regular expressions more readable. There is no special syntax for using them, they're just like single characters so here's a table.

<table>
<tr><th>Metacharacter</th><th>Represents</th></tr>
<tr><td><code>\t</code></td><td>a tab</td></tr>
<tr><td><code>\s</code></td><td>any space</td></tr>
<tr><td><code>\w</code></td><td>any word character (it is the same as <code>[aA-zZ0-9_]</code>)</td></tr>
<tr><td><code>\d</code></td><td>any digit (it is the same as <code>[0-9]</code>)</td></tr>
<tr><td><code>\W</code></td><td>any non-word char</td></tr>
<tr><td><code>\D</code></td><td>any non-digit char</td></tr>
</table>

<div class="keypoints" markdown="1">
#### Key Points
* Operators specify patterns that simplify regular expressions
* Operators are simply characters that specify other more general patterns
* Metacharacters are special pairs of characters that denote classes of characters
* Anchors match the pattern _only_ when it occurs at the start or end of a string
</div>

<div class="challenge" markdown="1">
#### Review

So now we have enough knowledge to try a quick quiz. What does this pattern match `(wo.+d)` return, when applied to this string `How much wood, would a woodchuck chuck?`, that is to say what does this print out:

~~~
m = re.search('(wo.+d)', "How much wood, would a woodchuck chuck?")
print m.group(1)
~~~
{:class="in"}

1.   `wood`
2.   `would, would a wood`
3.   `would, would`

</div>

