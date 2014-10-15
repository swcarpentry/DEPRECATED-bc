---
layout: lesson
root: ../..
title: One last Wrinkle
---

Let's have one last look at the function we wrote to extract data from
lab notebooks:

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

We can make it easier to add new patterns to this function by making it
more declarative. The trick is to combine the regular expressions and
the IDs of the groups we want to return:

~~~
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
~~~
{:class="in"}

Each entry in the list `patterns` has two parts: a regular expression,
and then the indices of the group that will contain the year, month,
day, site, and reading if that pattern matches. The loop tries the
regular expressions in `patterns` one by one. As soon as a pattern
matches it returns the matched groups, permuting them according to the
indices so that the data always comes back in the same order. To handle
the format in Notebook \#3, we just add one line to this table:

~~~
['([A-Z][a-z]+) ([0-9]{1,2}) ([0-9]{4}) \\((.+)\\) (.+)', 3, 1, 2, 4, 5]
~~~

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

### More Tools

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

~~~
Granger's work on graphs \cite{dd-gr2007,gr2009},
particularly ones obeying Snape's Inequality
\cite{ snape87 } (but see \cite{quirrell89}),
has opened up new lines of research.  However,
studies at Unseen University \cite{stibbons2002,
stibbons2008} highlight several dangers.
~~~

Citations in LaTeX are written using `\cite{…}`, with cross-reference
labels in the curly braces. A single citation can include two or more
labels separated by commas. There may be white space before or after
labels or line breaks where a citation is split across two lines, and
there can be multiple citations per line.

Our first idea is to use a group to capture everything inside the curly
braces following the word 'cite':

~~~
m = re.search('cite{(.+)}', 'a \\cite{X} b')
print m.groups()
~~~
{:class="in"}
~~~
('X',)
~~~
{:class="out"}

It seems to work in one simple case, but what if there are multiple
citations on a single line?

~~~
m = re.search('cite{(.+)}', 'a \\cite{X} b \\cite{Y} c')
print m.groups()
~~~
{:class="in"}
~~~
('X} b \\cite{Y',)
~~~
{:class="out"}

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

~~~
m = re.search('cite{([^}]+)}', 'a \\cite{X} b')
print m.groups()
~~~
{:class="in"}
~~~
('X,)
~~~
{:class="out"}

This works for a single citation: all we've done is change '.' to the
negated set. What about multiple citations on a single line?

~~~
m = re.search('cite{([^}]+)}', 'a \\cite{X} b \\cite{Y} c')
print m.groups()
~~~
{:class="in"}
~~~
('X,)
~~~
{:class="out"}

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

~~~
print re.findall('cite{([^}]+)}', 'a \\cite{X} b \\cite{Y} c')
~~~
{:class="in"}
~~~
['X', 'Y']
~~~
{:class="out"}

It seems to produce the right output—not bad for a 7-character change.
What about spaces in citations?

~~~
print re.search('cite{([^}]+)}', 'a \\cite{ X} b \\cite{Y } c').groups()
~~~
{:class="in"}
~~~
[' X', 'Y ']
~~~
{:class="out"}

The good news is, nothing breaks. The bad news is, the spaces are saved
by `findall`, which isn't really what we want. We could tidy this up
after the fact using `string.strip`, but let's modify the pattern
instead:

~~~
print re.findall('cite{\\s*([^}]+)\\s*}', 'a \\cite{ X} b \\cite{Y } c')
~~~
{:class="in"}
~~~
['X', 'Y ']
~~~
{:class="out"}

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

~~~
print re.findall('cite{\\s*\\b([^}]+)\\b\\s*}', 'a \\cite{ X} b \\cite{Y } c')
~~~
{:class="in"}
~~~
['X', 'Y']
~~~
{:class="out"}

It works! check this last example: in the PowerPoint, there's still a
space before the 'X' The change is to put `'\b'` after the first
unwanted spaces, and before the last ones. Since the curly braces around
the citation labels are also non-word characters, the pattern matches
even if there aren't any opening or trailing spaces.

The last hurdle is to handle multiple labels inside a single pair of
curly braces. The pattern we've built so far doesn't explode when there
are two or more labels, and it even handles spaces after the commas, but
it returns all those labels as a single lump of text:

~~~
print re.findall('cite{\\s*\\b([^}]+)\\b\\s*}', '\\cite{X,Y} ')
~~~
{:class="in"}
~~~
['X,Y']
~~~
{:class="out"}

~~~
print re.findall('cite{\\s*\\b([^}]+)\\b\\s*}', '\\cite{X, Y, Z} ')
~~~
{:class="in"}
~~~
['X, Y, Z']
~~~
{:class="out"}

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

~~~
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
~~~
{:class="in"}
~~~
set([])
~~~
{:class="out"}

Now let's write our function. For readability's sake, we'll put our
patterns at the top and give them memorable names. Inside the function,
we'll pull out all the citations using the first pattern, then split
each result everywhere there's a comma with optional spaces before or
after it. We'll stuff all the results into a set, and return that. If no
matches were found, that set will be empty.

~~~
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
~~~
{:class="in"}

We can use one more trick from the regular expression library to make
this function more efficient. Instead of turning the regular expression
into a finite state machine over and over again, we can compile the
regular expression and save the resulting object:

~~~
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
~~~
{:class="in"}

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

~~~
if __name__ == '__main__':
    test = '''\
Granger's work on graphs \cite{dd-gr2007,gr2009},
particularly ones obeying Snape's Inequality
\cite{ snape87 } (but see \cite{quirrell89}),
has opened up new lines of research.  However,
studies at Unseen University \cite{stibbons2002,
stibbons2008} highlight several dangers.'''

    print get_citations(test)
~~~
{:class="in"}
~~~
set(['gr2009', 'stibbons2002', 'dd-gr2007', 'stibbons2008',
     'snape87', 'quirrell89'])
~~~
{:class="out"}

Finally, if we are going to compile our regular expressions, we can make
them even easier to understand by using *verbose mode* to add comments.
Verbose mode tells Python to ignore whitespace and comments in the
regular expression, which lets us write patterns like this:

~~~
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
~~~
{:class="in"}

Documenting patterns like this makes them much easier to fix and extend.

