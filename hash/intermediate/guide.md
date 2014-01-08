---
layout: lesson
root: ../..
title: Instructor's Guide
level: intermediate
---
The ostensible goal of this set of lessons is
to introduce learners to non-linear data structures.
Most have ony ever seen arrays or lists,
i.e.,
things that are accessed using sequential numeric indices.
Sets and dictionaries are usually their first exposure to
accessing content by value rather than by location,
and to the bigger idea that
there are lots of other data structures they might want to learn about.

These lessons also introduce JSON as a general-purpose data format
that requires less effort to work with than flat text or CSV.
We discuss its shortcomings as well as its benefits
to help learners see what forces are at play
when designing and/or choosing data representations.

Finally,
these lessons are our first chance to introduce the idea of computational complexity
via back-of-the-envelope calculations of
how the number of steps required to look things up in an unordered list
grows with the number of things being looked up.
The discussion of hash tables is also good preparation for
discussion of relation databases,
but isn't required.

Teaching Notes
--------------

*   Start with sets:
    they're a familiar concept,
    there's no confusion between keys and values,
    and they are enough to motivate discussion of hash tables.

*   Python's requirement that entries in hash-based data structures be immutable
    only makes sense once the mechanics of hash tables are explained.
    Terms like "hash codes" and "hash function" also come up
    in error messages and Python's documentation,
    so learners are likely to become confused without some kind of hand-waving overview.
    Tuples are also easy to explain as
    "how to create immutable multi-part keys",
    and it's easy to explain why entries can't be looked up by parts
    (e.g., why a tuple containing a first and a last name can't be looked up by last name only)
    in terms of hash functions.

*   Finally, explaining why hash tables are fast is a good first encounter with
    the idea of "big oh" complexity.

*   Once sets have been mastered,
    dictionaries can be explained as "sets with extra information attached to each entry".
    The canonical example&mdash;counting things&mdash;shows
    why that "extra information" is useful.

*   Use the nanotechnology inventory example to re-emphasize how code is built top-down
    by writing code as if desired functions existed,
    then filling them in.
