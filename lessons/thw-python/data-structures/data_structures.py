
# <markdowncell>
# # Compound Data Types: Lists, Dictionaries, Sets, Tuples, and Reading Files
# 
# * * * * *
# <markdowncell>
# **Based on lecture materials by Milad Fatenejad, Joshua R. Smith, and Will Trimble**
# 
# Python would be a farily useless language if it weren't for the compound data types. The main two are **lists** and **dictionaries**, but I'll mention **sets** and **tuples** as well. I'll also go over reading text data from files. 
# 
# <markdowncell>
# ## Lists
# 
# A list is an ordered, indexable collection of data. Lets say you have collected some current and voltage data that looks like this:
# 
# ```
# voltage:
# -2.0
# -1.0
# 0.0
# 1.0
# 2.0
# 
# current:
# -1.0
# -0.5
# 0.0
# 0.5
# 1.0
# ```
# 
# So you could put that data into lists like
# 
# <codecell>
voltageList = [-2.0, -1.0, 0.0, 1.0, 2.0]
currentList = [-1.0, -0.5, 0.0, 0.5, 1.0]

# <markdowncell>
# obviously voltageList is of type list:
# 
# <codecell>
type(voltageList)
# <markdowncell>
# 
# Python lists have the charming (annoying?) feature that they are indexed from zero. Therefore, to find the value of the first item in voltageList:
# 
# <codecell>
voltageList[0]
# <markdowncell>
# 
# And to find the value of the third item

# <codecell>
voltageList[2]
# <markdowncell>

# Lists can be indexed from the back using a negative index. The last item of currentList
# 
# <codecell>
currentList[-1]
# <markdowncell>
# 
# and the next-to-last
# <codecell>
currentList[-2]
# <markdowncell>
# 
# You can "slice" items from within a list. Lets say we wanted the second through fourth items from voltageList
# <codecell>
voltageList[1:4]
# <markdowncell>
# 
# Or from the third item to the end
# <codecell>
voltageList[2:]
# <markdowncell>
# 
# and so on.
# 
# <markdowncell>
# ### Append and Extend
# 
# Just like strings have methods, lists do too.
# 
# <codecell>
dir(list)
# <markdowncell>
# 
# One useful method is append. Lets say we want to stick the following data on the end of both our lists.
# 
# ```
# voltage:
# 3.0
# 4.0
# 
# current:
# 1.5
# 2.0
# ```
# 
# If you want to append items to the end of a list, use the append method.
# 
# <codecell>
voltageList.append(3.)
# <codecell>
voltageList.append(4.)
# <codecell>
voltageList
# <markdowncell>
# 
# You can see how that approach might be tedious in certain cases. If you want to concatenate a list onto the end of another one, use extend.
# 
# <codecell>
currentList.extend([1.5, 2.0])
# <codecell>
currentList
# <markdowncell>
# ### Length of Lists
# 
# Sometimes you want to know how many items are in a list. Use the len command.
# 
# <codecell>
len(voltageList)
# <markdowncell>
# 
# ### Heterogeneous Data
# 
# Lists can contain hetergeneous data.
# 
# <codecell>
dataList = ["experiment: current vs. voltage", \
    "run", 47, \
    "temperature", 372.756, \
    "current", [-1.0, -0.5, 0.0, 0.5, 1.0], \
    "voltage", [-2.0, -1.0, 0.0, 1.0, 2.0]]
# <codecell>
print dataList
# <markdowncell>
 
# We've got strings, ints, floats, and even other lists in there. The slashes are there so we can continue on the next line. They aren't necessary but they can sometimes make things look better.
# 
# ## Assigning Variables to Other Variables
# 
# Something that might cause you headaches in the future is how python deals with assignment of one variable to another. When you set a variable equal to another, both variables point to the same thing. Changing the first one ends up changing the second. Be careful about this fact.
# 
# <codecell>
a = [1,2]
# <codecell>
b = a
# <codecell>
a.append(10)
# <codecell>
b
# <markdowncell>
# 
# There's a ton more to know about lists, but lets press on. Check out Dive Into Python or the help documentation for more info.
# 
# ## Reading From Files
# 
# At this point it is useful to take a detour regarding files. Lets say you have a file with some current and voltage data and some metadata.
# 
# ```
# data.dat:
# 
# experiment: current vs. voltage
# run: 47
# temperature: 372.756
# current: [-1.0, -0.5, 0.0, 0.5, 1.0]
# voltage: [-2.0, -1.0, 0.0, 1.0, 2.0]
# ```
# 
# We can read this data into a list type variable pretty easily.
# 
# <codecell>
f = open("data.dat")
# <codecell>
ivdata = f.readlines()
# <codecell>
f.close()
# <codecell>
ivdata
# <markdowncell>
# 
# Right now the data in ivdata isn't in a particularly useful format, but you can imagine that with some additional programming we could straighten it out. We will eventually do that.
# 
# ## Tuples
# 
# Tuples are another of python's basic compound data types that are almost like lists. The difference is that a tuple is immutable; once you set the data in it, the tuple cannot be changed. You define a tuple as follows.
# 
# <codecell>
tup = ("red", "white", "blue")
 
# <codecell>
type(tup)
# <markdowncell>
# 
# You can slice and index the tuple exactly like you would a list. Tuples are used in the inner workings of python, and a tuple can be used as a key in a dictionary, whereas a list cannot as we will see in a moment.
# 
# See if you can retrieve the third element of **tup**:
# <codecell>
 
# <markdowncell>
# ## Sets
# 
# Most introductary python courses do not go over sets this early (or at all), but I've found this data type to be useful. The python set type is similar to the idea of a mathematical set: it is an unordered collection of unique things. Consider:
# 
# <codecell>
fruit = set(["apple", "banana", "pear", "banana"]) 
# <markdowncell>
#You have to use a list to create a set.

# 
# Since sets contain only unique items, there's only one banana in the set fruit.
# 
# You can do things like intersections, unions, etc. on sets just like in math. Here's an example of an intersection of two sets (the common items in both sets).
# 
# <codecell>
firstBowl = set(["apple", "banana", "pear", "peach"])
 
# <codecell>
secondBowl = set(["peach", "watermelon", "orange", "apple"])
 
# <codecell>
set.intersection(firstBowl, secondBowl)
# <codecell>
set(['apple', 'peach'])

# <markdowncell>
# 
# You can check out more info using the help docs. We won't be returning to sets, but its good for you to know they exist.
# 
# ## Dictionaries
# 
# Recall our file data.dat which contained our current-voltage data and also some metadata. We were able to import the data as a list, but clearly the list type is not the optial choice for a data model. The dictionary is a much better choice. A python dictionary is a collection of key, value pairs. The key is a way to name the data, and the value is the data itself. Here's a way to create a dictionary that contains all the data in our data.dat file in a more sensible way than a list.
# 
# <codecell>
dataDict = {"experiment": "current vs. voltage", \
                    "run": 47, \
                    "temperature": 372.756, \
                    "current": [-1.0, -0.5, 0.0, 0.5, 1.0], \
                    "voltage": [-2.0, -1.0, 0.0, 1.0, 2.0]}
# <codecell>
print dataDict

# <markdowncell>
# 
# This model is clearly better because you no longer have to remember that the run number is in the second position of the list, you just refer directly to "run":
# 
# <codecell>
dataDict["run"]
# <markdowncell>
# 
# If you wanted the voltage data list:
# 
# <codecell>
dataDict["voltage"]
# <markdowncell>
# 
# Or perhaps you wanted the last element of the current data list
# 
# <codecell>
dataDict["current"][-1]
# <markdowncell>
# 
# Once a dictionary has been created, you can change the values of the data if you like.
# 
# <codecell>
dataDict["temperature"] = 3275.39
# <markdowncell>
# 
# You can also add new keys to the dictionary.  Note that dictionaries are indexed with square braces, just like lists--they look the same, even though they're very different.
# 
# <codecell>
dataDict["user"] = "Johann G. von Ulm"
# <markdowncell>
# 
# Dictionaries, like strings, lists, and all the rest, have built-in methods. Lets say you wanted all the keys from a particular dictionary.
# 
# <codecell>
dataDict.keys()
# <markdowncell>
# 
# also, values
# 
# <codecell>
dataDict.values()
# <markdowncell>
# 
# The help documentation has more information about what dictionaries can do.
# 
# Its worth mentioning that the value part of a dictionary can be any kind of data, even another dictionary, or some complex nested structure. The same is true about a list: they can contain complex data types.
# 
# Since tuples are immutable, they can be used as keys for dictionaries. Lists are mutable, and therefore cannot.
# 
# When you design software in python, most data will end up looking either like a list or a dictionary. These two data types are very important in python and you'll end up using them all the time.
