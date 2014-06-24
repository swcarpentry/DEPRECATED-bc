# Introduction to Regular Expressions

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