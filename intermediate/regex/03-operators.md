# Operators


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