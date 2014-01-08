---
layout: lesson
root: ../..
title: JSON
level: intermediate
---
<div class="objectives" markdown="1">
## Objectives

*   FIXME
</div>

## Lesson

The example above used two data file formats: one for storing molecular
formulas, the other for storing inventory. Both formats were specific to
this application, which means we needed to write, debug, document, and
maintain functions to handle them. Those functions weren't particularly
difficult to create, but they still took time to create, and if anyone
ever wants to read our files in Java, MATLAB, or Perl, they'll have to
write equivalent functions themselves.

A growing number of programs avoid these problems by using a flexible
data format called [JSON](glossary.html#json), which stands for
"JavaScript Object Notation". Despite the name, it is a
language-independent way to store nested data structures made up of
strings, numbers, Booleans, lists, dictionaries, and the special value
`null` (equivalent to Python's `None`)—in short, the basic data types
that almost every language supports. For example, let's convert a
dictionary of scientists' birthdays to a string:

~~~~ {src="setdict/json_first.py"}
>>> import json
>>> birthdays = {'Curie' : 1867, 'Hopper' : 1906, 'Franklin' : 1920}
>>> as_string = json.dumps(birthdays)
>>> print as_string
{"Curie": 1867, "Hopper": 1906, "Franklin": 1920}
>>> print type(as_string)
<type 'str'>
~~~~

`json.dumps` doesn't seem to do much, but that's kind of the point: the
textual representation of the data structure looks pretty much like what
a programmer would type in to re-create it. The advantage is that this
representation can be saved in a file:

~~~~ {src="setdict/json_second.py"}
>>> import json
>>> original = {'Curie' : 1867, 'Hopper' : 1906, 'Franklin' : 1920}
>>>
>>> writer = open('/tmp/example.json', 'w')
>>> json.dump(original, writer)
>>> writer.close()
>>>
>>> reader = open('/tmp/example.json', 'r')
>>> duplicate = json.load(reader)
>>> reader.close()
>>>
>>> print 'original:', original
original: {'Curie': 1867, 'Hopper': 1906, 'Franklin': 1920}
>>> print 'duplicate:', duplicate
duplicate: {u'Curie': 1867, u'Hopper': 1906, u'Franklin': 1920}
>>> print 'original == duplicate:', original == duplicate
original == duplicate: True
>>> print 'original is duplicate:', original is duplicate
original is duplicate: False
~~~~

As the example above shows, saving and loading data is as simple as
opening a file and then calling one function. The data file holds what
we'd type in to create the data in a program:

    $ cat /tmp/example.json
    {"Curie": 1867, "Hopper": 1906, "Franklin": 1920}

which makes it easy to edit by hand.

How is this different in practice from what we had? First, our inventory
file now looks like this:

~~~~ {src="setdict/inventory.json"}
{"He" : 1,
 "H" : 4,
 "O" : 3}
~~~~

while our formulas files look like:

~~~~ {src="setdict/formulas.json"}
{"helium"   : {"He" : 1},
 "water"    : {"H" : 2, "O" : 1},
 "hydrogen" : {"H" : 2}}
~~~~

Those aren't as intuitive for non-programmers as the original flat text
files, but they're not too bad. The worst thing is the lack of comments:
unfortunately—very unfortunately—the JSON format doesn't support them.
(And note that JSON requires us to use a double-quote for strings:
unlike Python, we cannot substitute single quotes.)

The good news is that given files like these, we can rewrite our program
as:

~~~~ {src="setdict/nanotech_json.py"}
'''Calculate how many molecules of each type can be made with the atoms on hand.'''

import sys
import json

def main(argv):
    '''Main driver for program.'''
    inventory = read_data(argv[1])
    formulas = read_data(argv[2])
    counts = calculate_counts(inventory, formulas)
    show_counts(counts)

def read_data(filename):
    '''Read a JSON-formatted data file.'''
    reader = open(filename, 'r')
    result = json.load(reader)
    reader.close()
    return result

def calculate_counts(inventory, formulas):
    ...as before...

def dict_divide(inventory, molecule):
    ...as before...

def show_counts(counts):
    ...as before...

if __name__ == '__main__':
    main(sys.argv)
~~~~

The two functions that read formula and inventory files have been
replaced with a single function that reads JSON. Nothing else has to
change, because the data structures loaded from the data files are
exactly what we had before. The end result is 51 lines long compared to
the 80 we started with, a reduction of more than a third.

### Nothing's Perfekt

JSON's greatest weakness isn't its lack of support for comments, but the
fact that it doesn't recognize and manage aliases. Instead, each
occurrence of an aliased structure is treated as something brand new
when data is being saved. For example:

    >>> inner = ['name']
    >>> outer = [inner, inner] # Creating an alias.
    >>> print outer
    [['name'], ['name']]
    >>> print outer[0] is outer[1]
    True
    >>> as_string = json.dumps(outer)
    >>> duplicate = json.loads(as_string)
    >>> print duplicate
    [[u'name'], [u'name']]
    >>> print duplicate[0] is duplicate[1]
    False

[Figure 9](#f:json_alias) shows the difference between the original data
structure (referred to by `outer`) and what winds up in `duplicate`. If
aliases might be present in our data, and it's important to preserve
their structure, we must either record the aliasing ourself (which is
tricky), or use some other format. Luckily, a lot of data either doesn't
contain aliases, or the aliasing in it isn't important.

![Aliasing in JSON](setdict/json_alias.png)

Figure 9: Aliasing in JSON

### Summary

-   The JSON data format can represent arbitrarily-nested lists and
    dictionaries containing strings, numbers, Booleans, and `none`.
-   Using JSON reduces the code we have to write ourselves and improves
    interoperability with other programming languages.

### Challenges

1.  A friend of yours says, "I understand why flat text files are not
    ideal, but wouldn't it be better to use comma-separated values (CSV)
    than JSON? It's easier to read, and more programs support it." What
    example could you show your friend to explain JSON's advantages?
2.  `json.dump` has an extra parameter called `sort_keys`; its default
    value is `False`, but if it is `True`, then all dictionaries are
    printed with keys in sorted order. Explain why this option *isn't*
    `True` by default, and how setting it to `True` can be useful in
    testing.
3.  If we really do need to add comments to JSON files, how can we do it
    without altering the format?
4.  The birdwatching data from [an earlier section](#s:aggregation) was
    stored like this:

        2010-07-03    05:38    loon
        2010-07-03    06:02    goose
        2010-07-03    06:07    loon
        2010-07-04    05:09    ostrich
        2010-07-04    05:29    loon

    How would you represent this as JSON? If you rewrite the
    `early_bird.py` program (that finds the earliest time each bird was
    seen) so that it uses your JSON format, how much code do you save?

Phylogenetic Trees
------------------

### Learning Objectives

-   That many "matrix" problems may be best solved using dictionaries.
-   Why the values in multi-part keys should be ordered.

30 minutes (but only for proficient learners—instructors should *not*
try to present this material to learners with weaker programming
backgrounds).

As Theodosius Dobzhansky said almost a century ago, nothing in biology
makes sense except in the light of evolution. Since mutations usually
occur one at a time, the more similarities there are between the DNA of
two species, the more recently they had a common ancestor. We can use
this idea to reconstruct the evolutionary tree for a group of organisms
using a hierarchical clustering algorithm.

We don't have to look at the natural world very hard to realize that
some organisms are more alike than others. For example, if we look at
the appearance, anatomy, and lifecycles of the seven fish shown in
[Figure 10](#f:species_pairs), we can see that three pairs are closely
related. But where does the seventh fit in? And how do the pairs relate
to each other?

  ---------------------------------------------------- ----------------------------------------------------
  ![Pairing Up Species](setdict/species_pairs_1.png)   ![Pairing Up Species](setdict/species_pairs_2.png)
  ---------------------------------------------------- ----------------------------------------------------

Figure 10: Pairing Up Species

The first step is to find the two species that are most similar, and
construct their plausible common ancestor. We then pair two more, and
two more, and start joining pairs to individuals, or pairs with other
pairs. Eventually, all the organisms are connected. We can redraw those
connections as a tree, using the heights of branches to show the number
of differences between the species we're joining up ([Figure
11](#f:species_tree)).

![Pairing Up Species](setdict/species_pairs_3.png)

Figure 11: Tree of Life

Let's turn this into an algorithm:

    S = {all organisms}
    while S != {}:
      a, b = two closest entries in U
      p = common parent of {a, b}
      S = S - {a, b}
      S = S + {p}

Initially, the set S contains all the species we're interested in. Each
time through the loop, we find the two that are closest, create their
common parent, remove the two we just paired up from the set, and insert
the newly-created parent. Since the set shrinks by one element each time
(two out, one in), we can be sure this algorithm eventually terminates.

But how do we calculate the distance between an inferred parent and
other species? One simple rule is to use the average distance between
that other species and the two species that were combined to create that
parent. Let's illustrate it by calculating a phylogenetic tree for
humans, vampires, werewolves, and mermaids. The distances between each
pair of species is shown in [Figure 12](#f:species_tree) and in the
table below. (We only show the lower triangle because it's symmetric.)

  -------------------------------------------------------------- -------------------------------------------------------------- -------------------------------------------------------------- --------------------------------------------------------------
  ![Distances Between Species](setdict/species_distance_1.png)   ![Distances Between Species](setdict/species_distance_2.png)   ![Distances Between Species](setdict/species_distance_3.png)   ![Distances Between Species](setdict/species_distance_4.png)
  -------------------------------------------------------------- -------------------------------------------------------------- -------------------------------------------------------------- --------------------------------------------------------------

Figure 12: Distances Between Species

  ---------- ------- --------- ---------- ---------
             human   vampire   werewolf   mermaid
  human                                    
  vampire    13                            
  werewolf   5       6                     
  mermaid    12      15        29          
  ---------- ------- --------- ---------- ---------

The closest entries—i.e., the pair with minimum distance—are human and
werewolf. We replace this with a common ancestor, which we will call HW,
then set the distance between it and each other species X to be (HX +
WX)/2, i.e., the average of the human-to-X and werewolf-to-X distances.
This gives us a new table:

  --------- ------ --------- ---------
            HW     vampire   mermaid
  HW                          
  vampire   9.5               
  mermaid   20.5   15         
  --------- ------ --------- ---------

Repeating this step, we combine HW with V:

  --------- ------- ---------
            HWV     mermaid
  HWV                
  mermaid   17.75    
  --------- ------- ---------

and finally HWV with M.

We illustrated our algorithm with a triangular matrix, but the order of
the rows and columns is arbitrary. The matrix is really just a lookup
table mapping species to distances, and as soon as we think of lookup
tables, we should think of dictionaries. The keys are species—either the
ones we started with, or the ones we created—and the values are the
distances between them, so our original table becomes:

    {
        ('human',   'mermaid')  : 12,
        ('human',   'vampire')  : 13,
        ('human',   'werewolf') :  5,
        ('mermaid', 'vampire')  : 15,
        ('mermaid', 'werewolf') : 29,
        ('vampire', 'werewolf') :  6
    }

There is one trick here. Whenever we have a distance, such as that
between mermaids and vampires, we have to decide whether to use the key
`('mermaid', 'vampire')` or `('vampire', 'mermaid')` (or to record the
value twice, once under each key).

Let's start by setting up our test case and then calling a top-level
function to process our data:

~~~~ {src="setdict/phylogen.py"}
if __name__ == '__main__':

    species = {'human', 'mermaid', 'werewolf', 'vampire'}

    scores = {
        ('human',   'mermaid')  : 12,
        ('human',   'vampire')  : 13,
        ('human',   'werewolf') :  5,
        ('mermaid', 'vampire')  : 15,
        ('mermaid', 'werewolf') : 29,
        ('vampire', 'werewolf') :  6
    }

    order = main(species, scores)
    print order
~~~~

In a real program, of course, the data would be read in from a file, and
the set of actual species' names would be generated from it, but this
will do for now.

Next, let's translate our algorithm into something that could be
runnable Python:

~~~~ {src="setdict/phylogen.py"}
def main(species, scores):
    result = []
    while len(species) > 1:
        left, right = find_min_pair(species, scores)
        result.append(make_pair(left, right))
        species -= {left, right}
        make_new_pairs(species, scores, left, right)
        species.add(make_name(left, right))
    return result
~~~~

This is almost a direct translation of our starting point; the only
significant difference is that we're keeping the set of "active" species
in a set, and the scores in a dictionary. As species are combined, we
remove their names from the set and add a made-up name for their parent.
We never actually remove scores from the table; once the name of a
species is out of the set `species`, we'll never try to look up anything
associated with it in `scores` again.

The next step is to write `find_min_pair` to find the lowest score
currently in the table:

    def find_min_pair(species, scores):
        min_pair = None
        min_val = None
        for left in species:
            for right in species:
                if left < right:
                    this_pair = make_pair(left, right)
                    if (min_val is None) or (scores[this_pair] < min_val):
                        min_pair = this_pair
                        min_val = scores[this_pair]
        return min_pair

This function loops over all possible combinations of species names, but
only actually *uses* the ones that pass our ordering test (i.e., the
ones for which the first species name comes before the second species
name). If this is the first score we've looked at, or if it's lower than
a previously-seen score, we record the pair of species and the
associated score. When we're done, we return the pair of species.

The function that makes new entries for the table is fairly
straightforward as well. It just loops over all the active species,
averages the distances between them and the two species being combined,
and puts a new score in the table:

    def make_new_pairs(species, scores, left, right):
        for current in species:
            left_score = scores[make_pair(current, left)]
            right_score = scores[make_pair(current, right)]
            new_score = (left_score + right_score) / 2.0
            scores[make_pair(current, make_name(left, right))] = new_score

Finally, the `make_pair` and `make_name` functions are simply:

    def make_pair(left, right):
        if left < right:
            return (left, right)
        else:
            return (right, left)

    def make_name(left, right):
        return '<%s, %s>' % make_pair(left, right)

Let's try running the program:

    $ python phylogen.py
    [('human', 'werewolf'), ('<human, werewolf>', 'vampire'), ('<<human, werewolf>, vampire>', 'mermaid')]

This shows that humans and werewolves were combined first, that their
pairing was then combined with vampires, and that mermaids were added to
the cluster last. We obviously should do a lot more testing, but so far,
we seem to be on the right track.

<div class="keypoints" markdown="1">
## Key Points

*   FIXME
</div>

<div class="challenges" markdown="1">
## Challenges

*   FIXME
</div>
