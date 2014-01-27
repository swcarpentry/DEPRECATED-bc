---
layout: lesson
root: ../..
title: Aggregation
level: intermediate
---
<div class="objectives" markdown="1">
## Objectives

*   Recognize problems that can be solved by aggregating values.
*   Use dictionaries to aggregate values.
*   Explain why actual data values should be used as initializers rather
    than "impossible" values.
</div>

## Lesson

To see how useful dictionaries can be, let's switch tracks and do some
birdwatching. We'll start by asking how early in the day we saw each
kind of bird? Our data consists of the date and time of the observation,
the bird's name, and an optional comment:

    2010-07-03    05:38    loon
    2010-07-03    06:02    goose
    2010-07-03    06:07    loon
    2010-07-04    05:09    ostrich   # hallucinating?
    2010-07-04    05:29    loon
         …           …        …

Rephrasing our problem, we want the minimum of all the times associated
with each bird name. If our data was stored in memory like this:

    loon = ['05:38', '06:07', '05:20', ...]

the solution would simply be `min(loon)`, and similarly for the other
birds. However, we have to work with the data we have, so let's start by
reading our data file and creating a list of tuples, each of which
contains a date, time, and bird name as strings:

~~~~ {src="setdict/early_bird.py"}
def read_observations(filename):
    '''Read data, returning [(date, time, bird)...].'''

    reader = open(filename, 'r')
    result = []

    for line in reader:
        fields = line.split('#')[0].strip().split()
        assert len(fields) == 3, 'Bad line "%s"' % line
        result.append(fields)

    return result
~~~~

This function follows the pattern we've seen many times before. We set
up by opening the input file and creating an empty list that we'll
append records to. We then process each line of the file in turn.
Splitting the line on the `'#'` character and taking the first part of
the result gets rid of any comment that might be present; stripping off
whitespace and then splitting breaks the remainder into fields.

To prevent trouble later on, we check that there actually are three
fields before going on. (An industrial-strength version of this function
would also check that the date and time were properly formatted, but
we'll skip that for now.) Once we've done our check, we append the
triple containing the date, time, and bird name to the list we're going
to return.

Here's the function that turns that list of tuples into a dictionary:

~~~~ {src="setdict/early_bird.py"}
def earliest_observation(data):
    '''How early did we see each bird?'''

    result = {}
    for (date, time, bird) in data:
        if bird not in result:
            result[bird] = time
        else:
            result[bird] = min(result[bird], time)

    return result
~~~~

Once again, the pattern should by now be familiar. We start by creating
an empty dictionary, then use a loop to inspect each tuple in turn. The
loop explodes the tuple into separate variables for the date, time and
bird. If the bird's name is not already a key in our dictionary, this
must be the first time we've seen it, so we store the time we saw it in
the dictionary. If the bird's name is already there, on the other hand,
we keep the minimum of the stored time and the new time. This is almost
exactly the same as our earlier counting example, but instead of either
storing 1 or adding 1 to the count so far, we're either storing the time
or taking the minimum of it and the least time so far.

Now, what if we want to find out which birds were seen on particular
days? Once again, we are [aggregating](glossary.html#aggregation)
values, i.e., combining many separate values to create one new one.
However, since we probably saw more than one kind of bird each day, that
"new value" needs to be a collection of some kind. We're only interested
in which birds we saw, so the right kind of collection is a set. Here's
our function:

~~~~ {src="setdict/birds_by_date.py"}
def birds_by_date(data):
    '''Which birds were seen on each day?'''

    result = {}
    for (date, time, bird) in data:
        if date not in result:
            result[date] = {bird}
        else:
            result[date].add(bird)

    return result
~~~~

Again, we start by creating an empty dictionary, and then process each
tuple in turn. Since we're recording birds by date, the keys in our
dictionary are dates rather than bird names. If the current date isn't
already a key in the dictionary, we create a set containing only this
bird, and store it in the dictionary with the date as the key.
Otherwise, we add this bird to the set associated with the date. (As
always, we don't need to check whether the bird is already in that set,
since the set will automatically eliminate any duplication.)

Let's watch this function in action for the first few records from our
data:

  Input                          Dictionary
  ------------------------------ --------------------------------------------------------------------------
  *start*                        `{}`
  `2010-07-03  05:38  loon`      `{'2010-07-03' : {'loon'}}`
  `2010-07-03  06:02  goose`     `{'2010-07-03' : {'goose', 'loon'}}`
  `2010-07-03  06:07  loon`      `{'2010-07-03' : {'goose', 'loon'}}`
  `2010-07-04  05:09  ostrich`   `{'2010-07-03' : {'goose', 'loon'}, '2010-07-04' : {'ostrich'}}`
  `2010-07-04  05:29  loon`      `{'2010-07-03' : {'goose', 'loon'}, '2010-07-04' : {'ostrich', 'loon'}}`

For our last example, we'll figure out which bird we saw least
frequently—or rather, which *birds*, since two or more may be tied for
the low score. Forgetting that values may not be unique is a common
mistake in data crunching, and often a hard one to track down.

One way to solve this problem is to do two passes over our data.
Assuming we have already built a dictionary with bird names as keys and
counts as values, we can then loop over the entries to find the one(s)
with the minimum count:

    for bird in counts:
        if counts[bird] < least:
            least = counts[bird]

and then use another loop to build the set of bird names that share that
minimum value:

    result = set()
    for bird in counts:
        if counts[bird] == least:
            result.add(bird)

There's a flaw in this code, though: we have not initialized `least`.
It's not safe to set it to an arbitrary "large" value like 1000 because
we might actually have seen some kind of bird more than 1000 times.

Another choice that *will* work is to initialize it to the largest
integer the computer can represent, which Python stores in `sys.maxint`.
While it's still possible that we'll have seen everything more often
than this (imagine, for example, that we were counting atoms instead of
counting birds), there's no way we could get those observations into our
program.

The best choice is to initialize `least` from the data itself. We'll
start by setting it to the special value `None`, and then check for this
value inside our loop. If `least` is `None`, it means we haven't
processed any data yet, so we'll assign whatever value we're looking at
to `least`. After that, all the less-than comparisons will do the right
thing. The resulting code looks like this:

    least = None
    for bird in counts:
        if least is None:
            least = counts[bird]
        elif counts[bird] < least:
            least = counts[bird]

This is another common design pattern: if we don't know what range of
values our data might take on, we initialize variables with a "don't
know" marker (such as `None`), then replace those markers with actual
values as quickly as possible.

Looping over the data three times—one to build the dictionary of names
and counts, a second time to find the minimum count, and a third time to
find all the birds we've seen that often—is unnecessarily slow. We can
actually solve the problem by reading the data just once. The idea is to
keep a dictionary whose keys are counts and whose values are sets of
birds. Every time we see a bird, we take it out of its current set and
move it to the set associated with the next larger count. If that set
doesn't exist yet, we simply create it. Here's a trace of how our
dictionary-of-sets would change as we read our sample data:

  Input                          Dictionary
  ------------------------------ -------------------------------------------------------
  *start*                        `{}`
  `2010-07-03  05:38  loon`      `{1 : {'loon'}}`
  `2010-07-03  06:02  goose`     `{1 : {'goose', 'loon'}}`
  `2010-07-03  06:07  loon`      `{1 : {'goose'}, 2 : {'loon'}}`
  `2010-07-04  05:09  ostrich`   `{1 : {'goose', 'ostrich'}, 2 : {'loon'}}`
  `2010-07-04  05:29  loon`      `{1 : {'goose', 'ostrich'}, 2 : set(), 3 : {'loon'}}`

The last step is the tricky one: when we read our final record, and move
`'loon'` from the set associated with 2 to a new set associated with 3,
there aren't any birds left that have been seen twice. We could remove
the empty set associated with the number 2, but it's simpler to leave it
there (and probably more efficient, too, since we're likely to promote
either `'goose'` or `'ostrich'` shortly thereafter).

But wait a second: how do we know which set a particular bird is
currently in? We can't look it up in our dictionary of sets, because the
sets holding birds' names are the values in that dictionary, and we can
only look up dictionary entries by key. The solution is to use another
dictionary storing the current count for each bird, just as we did
before:

~~~~ {src="setdict/least_common_bird.py"}
def least_common_birds(data):
    '''Which birds have been seen least frequently?'''

    count_by_bird = {}
    bird_by_count = {}
    for (date, time, bird) in data:

        # Record how often each bird has been seen.
        if bird in count_by_bird:
            old_count = count_by_bird[bird]
            count_by_bird[bird] = old_count + 1
        else:
            old_count = 0
            count_by_bird[bird] = 1

        # Update the birds stored by count.
        if count_by_bird[bird] == 1:
            if 1 in bird_by_count:
                bird_by_count[1].add(bird)
            else:
                bird_by_count[1] = {bird}
        else:
            bird_by_count[old_count].remove(bird)
            bird_by_count[old_count + 1].add(bird)

    # Return the set associated with the least key.
    # Note that there are (at least) two bugs in this block of code.
    keys = count_by_bird.keys()
    keys.sort()
    least = keys[0]
    return by_bird[least]
~~~~

As the comment says, there are at least two bugs in the final block of
code responsible for return the set of birds seen least frequently. We
will tackle them in the challenges.

<div class="keypoints" markdown="1">
## Key Points

*   Use dictionaries to count things.
*   Initialize values from actual data instead of trying to guess what
    values could "never" occur.
</div>

<div class="challenges" markdown="1">
## Challenges

1.  Draw a blob-and-arrow diagram of the two dictionaries in
    `least_common_birds` and all the data they refer to after the
    following seven lines of data have been processed:

        2013-06-23    05:31    sparrow
        2013-06-25    06:19    robin
        2013-07-03    06:21    robin
        2013-07-17    05:28    cardinal
        2013-07-19    05:28    robin
        2013-07-19    05:29    penguin
        2013-07-19    05:30    penguin

2.  The last four lines of the `least_common_birds` function in our
    final example can return the incorrect answer in at least two cases.
    Explain what these are and what wrong answer is returned in each,
    and then fix the function so that it does the right thing in both
    cases.
3.  We have frequently used the idiom:

        if key in data:
            data[key] = data[key] + 1
        else:
            data[key] = 1

    to either update the value associated with a key, or insert a value
    if the key isn't present. We can instead write this as:

        data[key] = data.get(key, 0) + 1

    Rewrite the `least_common_birds` function in our final example to
    use this idiom wherever possible, and explain why we *can't*
    simplify it even further by writing:

        data.get(key, 0) += 1

4.  Modify `least_common_birds` so that it returns a list of all the
    birds that have been seen, sorted from least common to most common.
    (Birds that appeared with equal frequency should be sorted
    alphabetically by name.)
5.  Write a function `dict_subtract` that "subtracts" one dictionary
    mapping names to numbers from another. For example:

        assert dict_subtract({'X' : 3},          {'X' : 2})          == {'X' : 1}
        assert dict_subtract({'X' : 3, 'Y' : 2}, {'X' : 5, 'Z' : 1}) == {'X' : -2, 'Y' : 2, 'Z' : -1}
</div>
