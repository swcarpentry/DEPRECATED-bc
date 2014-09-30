---
layout: lesson
root: ../..
title: Under the Hood
---

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
