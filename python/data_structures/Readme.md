# Compound Data Types: Lists, Dictionaries, Sets, Tuples, and Reading Files

* * * * *

**Based on lecture materials by Milad Fatenejad, Joshua R. Smith, and Will
Trimble**

Python would be a farily useless language if it weren't for the compound
data types. The main two are lists and dictionaries, but I'll mention sets
and tuples as well. I'll also go over reading text data from files. 

## Lists

A list is an ordered, indexable collection of data. Lets say you have
collected some current and voltage data that looks like this:

```
voltage:
-2.0
-1.0
0.0
1.0
2.0

current:
-1.0
-0.5
0.0
0.5
1.0
```

So you could put that data into lists like

```python
In [1]: voltageList = [-2.0, -1.0, 0.0, 1.0, 2.0]

In [2]: currentList = [-1.0, -0.5, 0.0, 0.5, 1.0]
```

obviously voltageList is of type list:

```python
In [3]: type(voltageList)
Out[3]: <type 'list'>
```

Python lists have the charming (annoying?) feature that they are indexed
from zero. Therefore, to find the value of the first item in voltageList:

```python
In [4]: voltageList[0]
Out[4]: -2.0
```

And to find the value of the third item

```python
In [5]: voltageList[2]
Out[5]: 0.0
```

Lists can be indexed from the back using a negative index. The last item of
currentList

```python
In [6]: currentList[-1]
Out[6]: 1.0
```

and the next-to-last

```python
In [7]: currentList[-2]
Out[7]: 0.5
```

You can "slice" items from within a list. Lets say we wanted the second
through fourth items from voltageList

```python
In [8]: voltageList[1:4]
Out[8]: [-1.0, 0.0, 1.0]
```

Or from the third item to the end

```python
In [9]: voltageList[2:]
Out[9]: [0.0, 1.0, 2.0]
```

and so on.

### Append and Extend

Just like strings have methods, lists do too.

```python
In [10] dir(list)
```

One useful method is append. Lets say we want to stick the following data
on the end of both our lists.

```
voltage:
3.0
4.0

current:
1.5
2.0
```

If you want to append items to the end of a list, use the append method.

```python
In [11]: voltageList.append(3.)

In [12]: voltageList.append(4.)

In [13]: voltageList
Out[13]: [-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0]
```

You can see how that approach might be tedious in certain cases. If you
want to concatenate a list onto the end of another one, use extend.

```python
In [14]: currentList.extend([1.5, 2.0])

In [15]: currentList
Out[15]: [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]
```

### Length of Lists

Sometimes you want to know how many items are in a list. Use the len command.

```python
In [16]: len(voltageList)
Out[16]: 7
```

### Heterogeneous Data

Lists can contain hetergeneous data.

```python
In [17]: dataList = ["experiment: current vs. voltage", \
   ....:             "run", 47, \
   ....:             "temperature", 372.756, \
   ....:             "current", [-1.0, -0.5, 0.0, 0.5, 1.0], \
   ....:             "voltage", [-2.0, -1.0, 0.0, 1.0, 2.0]]

```

We've got strings, ints, floats, and even other lists in there. The slashes
are there so we can continue on the next line. They aren't necessary but
they can sometimes make things look better.

## Assigning Variables to Other Variables

Something that might cause you headaches in the future is how python deals
with assignment of one variable to another. When you set a variable equal
to another, both variables point to the same thing. Changing the first one
ends up changing the second. Be careful about this fact.

```python
In [19]: a = [1,2]

In [20]: b = a

In [21]: a.append(10)

In [22]: b
Out[22]: [1, 2, 10]
```

There's a ton more to know about lists, but lets press on. Check out Dive
Into Python or the help documentation for more info.

## Reading From Files

At this point it is useful to take a detour regarding files. Lets say you
have a file with some current and voltage data and some metadata.

```
data.dat:

experiment: current vs. voltage
run: 47
temperature: 372.756
current: [-1.0, -0.5, 0.0, 0.5, 1.0]
voltage: [-2.0, -1.0, 0.0, 1.0, 2.0]
```

We can read this data into a list type variable pretty easily.

```python
In [1]: f = open("data.dat")

In [2]: ivdata = f.readlines()

In [3]: f.close()

In [4]: ivdata
Out[4]: 
['experiment: current vs. voltage\n',
 'run: 47\n',
 'temperature: 372.756\n',
 'current: [-1.0, -0.5, 0.0, 0.5, 1.0]\n',
 'voltage: [-2.0, -1.0, 0.0, 1.0, 2.0]\n',
 '\n']
```

Right now the data in ivdata isn't in a particularly useful format, but you
can imagine that with some additional programming we could straighten it
out. We will eventually do that.

## Tuples

Tuples are another of python's basic compound data types that are almost
like lists. The difference is that a tuple is immutable; once you set the
data in it, the tuple cannot be changed. You define a tuple as follows.

```python
In [1]: tup = ("red", "white", "blue")

In [2]: type(tup)
Out[2]: <type 'tuple'>
```

You can slice and index the tuple exactly like you would a list. Tuples are
used in the inner workings of python, and a tuple can be used as a key in a
dictionary, whereas a list cannot as we will see in a moment.

## Sets

Most introductary python courses do not go over sets this early (or at
all), but I've found this data type to be useful. The python set type is
similar to the idea of a mathematical set: it is an unordered collection of
unique things. Consider:

```python
In [3] fruit = set(["apple", "banana", "pear", "banana"]) #You have to use a list to create a set.
```

Since sets contain only unique items, there's only one banana in the set
fruit.

You can do things like intersections, unions, etc. on sets just like in
math. Here's an example of an intersection of two sets (the common items in
both sets).

```python
In [4]: firstBowl = set(["apple", "banana", "pear", "peach"])

In [5]: secondBowl = set(["peach", "watermelon", "orange", "apple"])

In [6]: set.intersection(firstBowl, secondBowl)
Out[6]: set(['apple', 'peach'])
```

You can check out more info using the help docs. We won't be returning to
sets, but its good for you to know they exist.

## Dictionaries

Recall our file data.dat which contained our current-voltage data and also
some metadata. We were able to import the data as a list, but clearly the
list type is not the optial choice for a data model. The dictionary is a
much better choice. A python dictionary is a collection of key, value
pairs. The key is a way to name the data, and the value is the data itself.
Here's a way to create a dictionary that contains all the data in our
data.dat file in a more sensible way than a list.

```python
In [7] dataDict = {"experiment": "current vs. voltage", \
                   "run": 47, \
                   "temperature": 372.756, \
                   "current": [-1.0, -0.5, 0.0, 0.5, 1.0], \
                   "voltage": [-2.0, -1.0, 0.0, 1.0, 2.0]}
```

This model is clearly better because you no longer have to remember that
the run number is in the second position of the list, you just refer
directly to "run":

```python
In [9]: dataDict["run"]
Out[9]: 47
```

If you wanted the voltage data list:

```python
In [10]: dataDict["voltage"]
Out[10]: [-2.0, -1.0, 0.0, 1.0, 2.0]
```

Or perhaps you wanted the last element of the current data list

```python
In [11]: dataDict["current"][-1]
Out[11]: 1.0
```

Once a dictionary has been created, you can change the values of the data
if you like.

```python
In [12]: dataDict["temperature"] = 3275.39
```

You can also add new keys to the dictionary.

```python
In [13]: dataDict["user"] = "Johann G. von Ulm"
```

Dictionaries, like strings, lists, and all the rest, have built-in methods.
Lets say you wanted all the keys from a particular dictionary.

```python
In [14]: dataDict.keys()
Out[14]: ['run', 'temperature', 'current', 'experiment', 'user', 'voltage']
```

also, values

```python
In [15]: dataDict.values()
Out[15]: 
[47,
 3275.39,
 [-1.0, -0.5, 0.0, 0.5, 1.0],
 'current vs. voltage',
 'Johann G. von Ulm',
 [-2.0, -1.0, 0.0, 1.0, 2.0]]
```

The help documentation has more information about what dictionaries can do.

Its worth mentioning that the value part of a dictionary can be any kind of
data, even another dictionary, or some complex nested structure. The same
is true about a list: they can contain complex data types.

Since tuples are immutable, they can be used as keys for dictionaries.
Lists are mutable, and therefore cannot.

When you architect software in python, most data will end up looking either
like a list or a dictionary. These two data types are very important in
python and you'll end up using them all the time.
