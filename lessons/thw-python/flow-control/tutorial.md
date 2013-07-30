---
layout: lesson
root: ../../..
title: Flow Control in Python
---
**Based on materials by Milad Fatenejad and Katy Huff**

Pasting into IPython
====================

This part of the lesson includes a lot of text, but it will be useful to
run it yourself in IPython.

To paste text from another application (i.e. these lecture notes) into
IPython :

1.  select text from the wiki
2.  copy with **ctrl+c**
3.  in IPython, type `%paste`

The code should paste and execute in IPython.

If you also type %autocall to turn autocall OFF, you may be able to
paste with **ctrl+v** though this won't work with all IPython builds.

Conditionals
============

A conditional (if statement) is some statement that in general says :
"When some boolean is true, do the following. Elsewise, do this other
thing."

Many equivalence test statements exist in Python that are similar in
other languages:

```python
i=1
j=2
i==j # i is equal to j : FALSE
i<j  # i is less than j
i<=j # i is less than or equal to j : TRUE
i>j  # i is greater than j
i>=j # i is greater than or equal to j : FALSE
i!=j # i is not equal to j : TRUE
```

However, python has other equivalence test statements that are fairly
unique to python. To check whether an object is contained in a list :

```python
beatle="John"
beatles=["George", "Ringo","John", "Paul"]
print beatle in beatles # is John one of the beatles? : TRUE
print "Katy" not in beatles # this is also TRUE. 
```

The `is` keyword tests if two variables refer to the same object. For example:

```python
a = 1234
b = 1234
a == b # True, they have the same value
a is b # False, are different objects
```

It is a common mistake to use `is` to test for equality between two objects, see the code below.  This only works for a small range of integers and strings in CPython, and is a side effect of the implementation that should **not** be relied upon.

```python
a = 1
b = 1
a is b # True - special case for 1
```

A correct use of `is` would be to compare objects like lists, for example the same list could be inserted into two different dictionaries. A comparison with `is` would reveal this:

```python
number_list = [1,2,4,8]
dict1 = {"thing_widths":number_list}
dict2 = {"item_costs":number_list}
dict1["thing_widths"] is dict2["item_costs"]  # True - this is the same list
```

Note that since the two dictionary values are actually the same object, modifying the list in one of the dictionaries will mean that the values in both dictionaries will change:

```python
print dict1, dict2
dict1["thing_widths"][0] = 222
print dict1, dict2
```

Conditionals (`if` statements) are also really easy to use in python. Take
a look at the following example:

```python
i = 4
sign = "zero"
if i < 0:
  sign = "negative"
elif i > 0:
  sign = "positive"
else:
  print "Sign must be zero"
  print "Have a nice day"
print sign
```

The behavior of this code snippet should be pretty clear, but there is
something peculiar. How does Python know where the if-statement ends?
Other languages, like FORTRAN, MatLab, and C/C++ all have some way of
delimiting blocks of code. For example, in MatLab you begin an if
statement with the word `if` and you end it with `end if`. In C/C++ you
delimit blocks with curly braces. Python uses **indentation** to delimit
code blocks. The **indentation** above is NOT just to make things look
pretty - it tells Python what the body of the `if`-statement is. This is
true when ever we create any code blocks, such as the bodies of loops,
functions or classes.

**Aside: Compact if-statement:**

Python has an easy to use `if`-syntax for setting the value of a variable.
Try entering this into IPython:

```python
i = 5
sign = "positive" if i > 0 else "negative"
```

While Loops
===========

Lets start by looking at while loops since they function like while
loops in many other languages. The example below takes a list of
integers and computes the product of each number in the list up to the
-1 element.

A while loop will repeat the instructions within itself until the
conditional that defines it is no longer true.

```python
mult = 1
sequence = [1, 5, 7, 9, 3, -1, 5, 3]
while sequence[0] != -1:
  mult = mult * sequence[0]
  del sequence[0]

print mult
```

Some new syntax has been introduced in this example.

-   On line 4, we compute the product of the elements just to make this
    more interesting.

-   On line 5, we use the `del` keyword to remove the first element of
    the list, shifting every element down one.

**Watch Out**

Since a while loop will continue until its conditional is no longer
true, a **poorly formed** while loop might repeat forever. For example :

```python
i=1
print "Well, there's egg and bacon, egg and spam, egg bacon and"
while i == 1:
  print "spam "
print "or Lobster Thermidor a Crevette with a mornay sauce served in a Provencale manner with shallots..." 
```

Since the variable `i` never changes within the while loop, we can
expect that the conditional, `i==1` will remain true forever and the
while loop will just go round and round, as if this restaurant offered
nothing but spam. (If you try this at home, please nory, or other iterable). However, sometimes you may need the index value at the same time, for example for some calculation. The `enumerate` function generates the integer index for you, which can be used instead of the `range` function. The following two loops are equivalent:

```python
data_list = [23,45,67]

for i in range(len(data_list)):
    print data_list[i], ' is item number ', i, ' in the list'

for i,d in enumerate(data_list):
    print d, ' is item number ', i, ' in the list'
```

###zip###

For iterating through multiple sequences, `zip` can be used to group them together to simultaneous pass through each sequence:

```python
run_numbers = [1,2,3,4]
run_times = [12.1, 33.0, 15.1, 22.9]
directions = ['North', 'South', 'East', 'NorthEast']

for n, t, d in zip(run_numbers, run_times, directions):
    print n, t, d
```

Final Example
=============

We've seen a lot so far. Lets work through a slightly lengthier example
together. I'll use some of the concepts we already saw and introduce a
few new concepts. To run the example, you'll need to locate a short file
containing phone numbers. The file can be found in your 
repository within the phonenums directory and is called phonenums.txt. 
Now we have to move IPython to that directory so it can find the
phonenums.txt file. You navigate within IPython in the same way that you
navigate in the shell, by entering "cd [path]" .

This example opens a text file containing a list of phone numbers. The
phone numbers are in the format \#\#\#-\#\#\#-\#\#\#\#, one to a line.
The example code loops through each line in the file and counts the
number of times each area code appears. The answer is stored in a
dictionary, where the area code is the key and the number of times it
occurs is the value.

```python

areacodes = {} # Create an empty dictionary
f = open("phonenums.txt") # Open the text file
for line in f: # iterate through the text file, one line at a time (think of the file as a list of lines)
    ac = line.split('-')[0] # Split phone number, first element is the area code
    if not ac in areacodes: # Check if it is already in the dictionary
        areacodes[ac] = 1 # If not, add it to the dictionary
    else:
        areacodes[ac] += 1 # Add one to the dictionary entry
  
print areacodes # Print the answer
```

Example : Iteritems
-------------------

Use the iteritems dictionary method in combination with a for loop to
print the keys/values of the areacodes dictionary one to a line. In
other words, the goal is to write a loop that prints:

    203 4
    800 4
    608 8
    773 3

This example is a little tricky to figure out, but give it a shot.
