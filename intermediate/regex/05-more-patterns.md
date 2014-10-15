---
layout: lesson
root: ../..
title: More Patterns
---

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
sequence](../../gloss.html#escape-sequence): just as we use the
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

~~~
# find '()' in text
m = re.search('\\(\\)', text)
    ⋮    ⋮    ⋮
~~~

The file has two backslashes, an open parenthesis, two backslashes, and
a close parenthesis inside quotes. When Python reads that file in, it turns the two-character sequence
'\\\\' into a single '\\' character in the string in memory. This is the
first round of escaping. When we hand that string '\\(\\)' to the regular expression library, it
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

### Regular Expressions and Newlines

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
