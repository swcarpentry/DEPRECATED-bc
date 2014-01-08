---
layout: lesson
root: ../..
title: Nanotech Inventory
level: intermediate
---
<div class="objectives" markdown="1">
## Objectives

*   Parse flat text files to create dictionaries.
*   Explain the pros and cons of JSON compared to flat text.
*   Create and manipulate nested dictionaries.
*   Explain the similarities and differences between nested dictionaries
    and nested lists.
</div>

## Lesson

We can now solve Fan's original nanotech inventory problem. As explained
in the introduction, our goal is to find out how many molecules of
various kinds we can make using the atoms in our warehouse. The number
of molecules of any particular type we can make is limited by the
scarcest atom that molecule requires. For example, if we have five
nitrogen atoms and ten hydrogen atoms, we can only make three ammonia
molecules, because we need three hydrogen atoms for each.

The formulas for the molecules we know how to make are stored in a file
like this:

    # Molecular formula file

    helium : He 1
    water : H 2 O 1
    hydrogen : H 2

and our inventory is stored in a file like this:

    # Atom inventory file

    He 1
    H 4
    O 3

Let's start by reading in our inventory. It consists of pairs of strings
and numbers, which by now should suggest using a dictionary for storage.
The keys will be atomic symbols, and the values will be the number of
atoms of that kind we currently have ([Figure
6](#f:nanotech_inventory)). If an atom isn't listed in our inventory,
we'll assume that we don't have any.

![Nanotech Inventory](setdict/nanotech_inventory.png)

Figure 6: Nanotech Inventory

What about the formulas for the molecules we know how to make? Once
again, we want to use strings—the names of molecules—as indices, which
suggests a dictionary. Each of its values will be something storing
atomic symbols and the number of atoms of that type in the molecule—the
same structure, in fact, that we're using for our inventory. [Figure
7](#f:nanotech_formulas) shows what this looks like in memory if the
only molecules we know how to make are water and ammonia.

![Storing Formulas](setdict/nanotech_formulas.png)

Figure 7: Storing Formulas

Finally, we'll store the results of our calculation in yet another
dictionary, this one mapping the names of molecules to how many
molecules of that kind we can make ([Figure 8](#f:nanotech_results)).

![Nanotech Results](setdict/nanotech_results.png)

Figure 8: Nanotech Results

The main body of the program is straightforward: it reads in the input
files, does our calculation, and prints the result:

~~~~ {src="setdict/nanotech.py"}
'''Calculate how many molecules of each type can be made with the atoms on hand.'''

import sys

if __name__ == '__main__':
    inventory = read_inventory(sys.argv[1])
    formulas = read_formulas(sys.argv[2])
    counts = calculate_counts(inventory, formulas)
    show_counts(counts)
~~~~

Reading the inventory file is simple. We take each interesting line in
the file, split it to get an atomic symbol and a count, and store them
together in a dictionary:

~~~~ {src="setdict/nanotech.py"}
def read_inventory(filename):
    '''Read inventory of available atoms.'''

    result = {}
    for line in read_lines(filename):
        name, count = line.split(' ')
        result[name] = int(count)

    return result
~~~~

For clarity's sake, we have used a helper function called `read_lines`
to remove blank lines and comments from our input:

~~~~ {src="setdict/nanotech.py"}
def read_lines(filename):
    '''Read lines from file, stripping out blank lines and comments.'''

    reader = open(filename, 'r')
    lines = []
    for line in reader:
        line = line.split('#')[0].strip()
        if line:
            lines.append(line)
    reader.close()

    return lines
~~~~

Using that same function, the function that reads in a file of molecular
formulas is only slightly more complex than the one that reads in
inventory:

~~~~ {src="setdict/nanotech.py"}
def read_formulas(filename):
    '''Read molecular formulas from file.'''

    result = {}                                        # 1
    for line in read_lines(filename):

        name, atoms = line.split(':')                  # 2
        name = name.strip()

        atoms = atoms.strip().split(' ')               # 3
        formula = {}
        for i in range(0, len(atoms), 2):              # 4
            formula[atoms[i]] = int(atoms[i+1])        # 5

        result[name] = formula                         # 6

    return result                                      # 7
~~~~

We start by creating a dictionary to hold our results (\#1). We then
split each interesting line in the data file on the colon ':' to
separate the molecule's name (which may contain spaces) from its formula
(\#2). We then split the formulas into a list of strings (\#3). These
alternate between atomic symbols and numbers, so in the inner loop
(\#4), we move forward through those values two elements at a time,
storing the atomic symbol and count in a dictionary (\#5). Once we're
done, we store that dictionary as the value for the molecule name in the
main dictionary (\#6). When we've processed all the lines, we return the
final result (\#7).

Now that we have all our data, it's time to calculate how many molecules
of each kind we can make. `inventory` maps atomic symbols to counts, and
so does `formulas[name]`, so let's loop over all the molecules we know
how to make and "divide" the inventory by each one:

~~~~ {src="setdict/nanotech.py"}
def calculate_counts(inventory, formulas):
    '''Calculate how many of each molecule can be made with inventory.'''

    counts = {}
    for name in formulas:
        counts[name] = dict_divide(inventory, formulas[name])

    return counts
~~~~

We say we're "dividing" the inventory by each molecule because we're
trying to find out how many of that molecule we can make without
requiring more of any particular atom than we actually have. (By
analogy, when we divide 11 by 3, we're trying to find out how many 3's
we can make from 11.) The function that does the division is:

~~~~ {src="setdict/nanotech.py"}
def dict_divide(inventory, molecule):
    '''Calculate how much of a single molecule can be made with inventory.'''

    number = None
    for atom in molecule:
        required = molecule[atom]
        available = inventory.get(atom, 0)
        limit = available / required
        if (number is None) or (limit < number):
            number = limit

    return number
~~~~

This function loops over all the atoms in the molecule we're trying to
build, see what limit the available inventory puts on us, and return the
minimum of all those results. This function uses a few patterns that
come up frequently in many kinds of programs:

1.  The first pattern is to initialize the value we're going to return
    to `None`, then test for that value inside the loop to make sure we
    re-set it to a legal value the first time we have real data. In this
    case, we could just as easily use -1 or some other impossible value
    as an "uninitialized" flag for `number`.
2.  Since we're looping over the keys of `molecule`, we know that we can
    get the value stored in `molecule[atom]`. However, that atom might
    not be a key in `inventory`, so we use `inventory.get(atom, 0)` to
    get either the stored value or a sensible default. In this case
    zero, the sensible default is 0, because if the atom's symbol isn't
    in the dictionary, we don't have any of it. This is our second
    pattern.
3.  The third is using calculate, test, and store to find a single
    value—in this case, the minimum—from a set of calculated values. We
    could calculate the list of available over required values, then
    find the minimum of the list, but doing the minimum test as we go
    along saves us having to store the list of intermediate values. It's
    probably not a noticeable time saving in this case, but it would be
    with larger data sets.

The last step in building our program is to show how many molecules of
each kind we can make. We could just loop over our result dictionary,
printing each molecule's name and the number of times we could make it,
but let's put the results in alphabetical order to make it easier to
read:

~~~~ {src="setdict/nanotech.py"}
def show_counts(counts):
    '''Show how many of each kind of molecule we can make.'''

    names = counts.keys()
    names.sort()
    for name in names:
        print name, counts[name]
~~~~

It's time to test our code. Let's start by using an empty inventory and
a single formula:

Inventory

Formulas

Output

    # inventory-00.txt

    # formulas-01.txt
    helium : He 1

There's no output, which is what we expect. Let's add one atom of helium
to our inventory:

Inventory

Formulas

Output

    # inventory-01.txt
    He 1

    # formulas-01.txt
    helium : He 1

    helium 1

That seems right as well. Let's add some hydrogen, but don't give the
program any formulas that use hydrogen:

Inventory

Formulas

Output

    # inventory-02.txt
    He 1
    H 4

    # formulas-01.txt
    helium : He 1

    helium 1

The output doesn't change, which is correct. Let's try adding the
formula for water, which does use hydrogen, but not providing any
oxygen:

Inventory

Formulas

Output

    # inventory-02.txt
    He 1
    H 4

    # formulas-02.txt
    helium : He 1
    water: H 2 O 1

    helium 1

As we hoped, there's no water in the output, but helium is still
appearing as it should. Let's add the formula for molecular hydrogen:

Inventory

Formulas

Output

    # inventory-02.txt
    He 1
    H 4

    # formulas-03.txt
    helium : He 1
    water: H 2 O 1
    hydrogen: H 2

    helium 1
    hydrogen 2

Sure enough, we can make two molecules of hydrogen (each of which uses
two atoms). Finally, let's put some oxygen in the warehouse:

Inventory

Formulas

Output

    # inventory-03.txt
    He 1
    H 4
    O 3

    # formulas-03.txt
    helium : He 1
    water: H 2 O 1
    hydrogen: H 2

    helium 1
    hydrogen 2
    water 2

That's right too: we can make two water molecules (because we don't have
enough hydrogen to pair with our three oxygen atoms). There are quite a
few other interesting tests still to run, but before we do that, we
should take another look at how we're storing our data on disk.

<div class="keypoints" markdown="1">
## Key Points

*   Whenever names are used to label things, consider using dictionaries
    to store them.
*   Use nested dictionaries to store hierarchical values (like molecule
    names and atomic counts).
</div>

<div class="challenges" markdown="1">
## Challenges

1.  Trace the behavior of `read_formulas` by showing the value of each
    variable each time line \#6 finishes executing when given the data
    file:

        helium  : He 1
        ammonia : N 1 H 3
        cyanide : H 1 C 1 N 1

                            `result`   `line`   `name`   `atoms`   `formula`
      --------------------- ---------- -------- -------- --------- -----------
      1) after "helium":                                           
      2) after "ammonia":                                          
      3) after "cyanide":                                          

2.  Can one dictionary be used as a key in another? I.e., is it possible
    to create the structure:

        { {'site' : 3, 'affinity' : 6} : 'sampled'}

    If so, give an example showing when this would be useful. If not,
    explain why not.

3.  A geographic information system stores the distance between survey
    points in a dictionary of dictionaries like this:

        dist = {
            'Left Bend'  : {'Sump Creek' : 25.6,
                            'Brents Bay' : 31.1,
                            'Ogalla'     :  4.0},
            'Sump Creek' : {'Brents Bay' : 17.5,
                            'Ogalla'     : 19.2},
            'Brents Bay' : {'Ogalla'     : 20.1}
        }

    Given this structure, what is the simplest Python function that will
    return the distance between any two survey points?

4.  Fan has inherited an activity log for an experimental project
    formatted as shown below:

        2012-11-30: Re-setting equipment.
        2012-12-12: First run with acidic reagants.
        2012-12-12: Re-ran acidic reagants.
        2012-12-14: Tried neutral reagants again.
        2013-02-05: Back to this stuff after writing up the CSRTI paper.
        2013-02-06: Trying basic reagants this time.

    He has written a function to translate this into a dictionary of
    dictionaries of sets, where the outer dictionary's keys are years
    (as strings), the inner dictionary's keys are months (also as
    strings), and the innermost sets are the days (strings again) for
    which there are comments. For example, the output of this function
    for the data sample above is supposed to be:

        {
            '2012' : {
                '11' : {'30'},
                '12' : {'12', '14'}
            },
            '2013' : {
                '02' : {'05', '06'}
            }
        }

    His function is:

        def extract_dates(filename):
            reader = open(filename, 'r')
            result = {}
            for line in reader:
                year, month, day = line.strip().split(' ', 1)[0].split('-')
                if year not in result:
                    pass # fill in 1
                if month not in result[year]:
                    pass # fill in 2
                pass # fill in 3
            reader.close()
            return result

    Fill in the three missing lines with a single statement each so that
    this function returns the right answer.
</div>
