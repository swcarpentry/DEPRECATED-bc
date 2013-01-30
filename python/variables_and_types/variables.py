# -*- coding: utf-8 -*-
# <nbformat>2.0</nbformat>

# <markdowncell>

# ## Python, iPython, and the basics
# 
# * * * * *
# 
# **Based on Lecture Materials By: Milad Fatenejad, Katy Huff, Tommy Guy, Joshua R. Smith, Will Trimble, and Many More**
# 
# <markdowncell>
# ## Introduction
# This lecture is on basic programming in python. In order to do the examples, we are going to use an environment called iPython notebook.  I expect this lecture to be interactive, so stop me at any point if you have questions. The correct power dynamic is that people are the masters and the machines are servants. The computer is a hammer; it exists to help us get things done.  We can hammer nails with the handle, with the claw of the hammer; some of us even hammer nails with bricks.  But when you learn what part of the hammer works best with nails, and have some experience swinging it, you spend less time worrying about the hammering and more time worrying about your furniture.
# 
# <markdowncell>
# So now would be a good time to roll out [PEP 20, The Zen of Python](http://www.python.org/dev/peps/pep-0020/)
# > Beautiful is better than ugly.  
# > Explicit is better than implicit.  
# > Simple is better than complex.  
# > Complex is better than complicated.  
# > Flat is better than nested.  
# > Sparse is better than dense.  
# > Readability counts.  
# > Special cases aren't special enough to break the rules.  
# > Although practicality beats purity.  
# > Errors should never pass silently.  
# > Unless explicitly silenced.  
# > In the face of ambiguity, refuse the temptation to guess.  
# > There should be one-- and preferably only one --obvious way to do it.  
# > Although that way may not be obvious at first unless you're Dutch.  
# > Now is better than never.   
# > Although never is often better than *right* now.  
# > If the implementation is hard to explain, it's a bad idea.  
# > If the implementation is easy to explain, it may be a good idea.  
# > Namespaces are one honking great idea -- let's do more of those!  
# 
# <markdowncell>
# Here is the reference material.
# 
# * [Dive into Python](http://www.diveintopython.net/toc/index.html)
# * [Software Carpentry's Python Lectures](http://software-carpentry.org/4_0/python/)
# * [IPython: A System for Interactive Scientific Computing](http://dx.doi.org/10.1109/MCSE.2007.53)
# * [How to Think Like a Computer Scientist](http://www.greenteapress.com/thinkpython/thinkpython.html)
# 
# <markdowncell>
# ## Lesson 1
# * print statements
# * variables
# * integers
# * floats
# * strings
# * types
# * type coersion
# * basic operations: add numbers, concatenate strings, basic data type functionality
# <markdowncell>
# ## Lesson 2
# * list
# * dictionary
# * set
# * tuple
# * file reading
# <markdowncell>
# ## Lesson 3
# * for loop
# * conditional (if) statements
# * while loops
# * iteration
# * writing to files
# <markdowncell>
# ## Lesson 4
# * methods
# * modules
# 
# 
# ## Python environments
# You can run python commands in a handful of ways; you can create executable scripts, you can run the python interpreter, you can run iPython, or you can run iPython notebook.  iPython is an alternative to the built-in Python interpreter with some nice features.  iPython notebook gives you interactive access to the python interpreter from within a browser window, and it allows you to save your commands as a "notebook".
# Let's give the built-in interpreter a spin just this once.  Open a **Terminal** window, which starts your default shell.  Type 

# <markdowncell>
# ``
# python 
# ``
# <markdowncell>
# And you should see python start up.  Type 
# <markdowncell>
# ``
# print "Fresh out of parrots"
# ``
# <markdowncell>
# Note the black-and-white wallpaper.
# Escape from python with 
# `` 
# quit()
# ``
# <markdowncell>
# ***
# 
# iPython has more useful features for interactive use than the standard python interpreter, but it works in the same way--it is interacitve, and you get there from the command line.  iPython notebook uses javascript to allow you to enter python commands in your browser and show you the result in the browser.  We'll use it from here on out.

# <codecell>

print "hello world"

# <markdowncell>

# ## Navigating in ipython notebook
# The box above is called the input cell; commands you put here will be fed to the python interpreter one at a time when you press **Shift-ENTER**.  
# The output of the command, or the error message, appears below the line you entered it on.
# The panel which may appear on the left has some notebook options; you can minimize the panel by double-clicking on the bar separating the windows. 
# <codecell>
print "Try and tell that to the young people"
print "of today--they won't believe you."
# <markdowncell>

# If you hit **ENTER** only, ipython notebook gives you another line in the current cell.  
# This allows you to compose multi-line commands and submit them to python all at once.  

# <markdowncell>

# Up and down arrows will allow you to move the cursor to different cells in the notebook, including these cells containing text (which you can edit in your browser).  
# Only the cells for which you press Shift-ENTER or Control-ENTER will be executed by the python interpreter.   

# <markdowncell>

# You can enter the same line over and over again into the interpreter.  It's wierd, but it's life. 

# <codecell>

i = 0

# <markdowncell>

# **Shift-ENTER** executes and moves to the next cell.  
# **Control-ENTER** executes and does *not* move to the next cell.  
# Try entering this cell a few times:  
# <codecell>

i = i + 1
print i

# <markdowncell>

# If you want to create new empty cells, it's three keys: **Shift-Enter**, **Control-M**, and then  **a**  This will insert more cells in the middle of the notebook. 

# <markdowncell>

# ## Getting Help
# 
# iPython has some nice help features. Let's say we want to know more about the integer data type. There are at least two ways to do this task:

# <codecell>

help(int)

# <markdowncell>

# which displays a scrolling text help, or

# <codecell>

int?

# <markdowncell>

# Which displays a shorter help summary in the magic pane at the bottom of the screen.  You can minimize the magic pane when it gets in your way.

# <markdowncell>

# If you wanted to see all the built-in commands available for something, use the *dir* command. Check out all of the methods of the object "Hello world", which are shared by all objects of the str type.

# <codecell>

dir("Hello world")

# <markdowncell>

# There's a method that looks important -- swapcase.  Let's see what it does:  

# <codecell>

"Hello world".swapcase()

# <markdowncell>

# Hrm.  Ahem.
# ## Executing code in files
# 
# If your code is in a file, you can execute it from the iPython shell with the **%run** command. Execute hello.py like so

# <codecell>

%run hellp.py

# <markdowncell>

# *Ooops.*  We misspelled **hello.py**, and python is giving us an error message.  Change the line above to hello.py, hit **Shift-ENTER**, and see what it does.

# <markdowncell>

# ## Clearing iPython
# 
# To clear everything from iPython, use the %reset command.

# <codecell>

mystring = "And three shall be the count." 
print mystring

# <codecell>

%reset

# <codecell>

print mystring

# <markdowncell>

# Note that the error message contains a recap of the input that caused the error (with an arrow, no less!)   It is objecting that **mystring** is not defined, since we just reset it.
# 
# ## Variables
# 
# All programming languages have variables, and python is no different. To create a variable, just name it and set it with the equals sign. One important caveat: variable names can only contain letters, numbers, and the underscore character. Let's set a variable.

# <codecell>

experiment = "current vs. voltage"

# <codecell>

print experiment

# <codecell>

voltage = 2

# <codecell>

current = 0.5

# <codecell>

print voltage, current

# <markdowncell>

# ## Types and Dynamic Typing
# 
# Like most programming languages, things in python are typed. The type refers to the type of data. We've already defined three different types of data in experiment, voltage, and current. The types are string, integer, and float. You can inspect the type of a variable by using the type command.

# <codecell>

type(experiment)

# <codecell>

type(voltage)

# <codecell>

type(current)

# <markdowncell>

# Python is a dynamically typed language (unlike, say, C++). If you know what that means, you may be feeling some fear and loathing right now. If you don't know what dynamic typing means, the next stuff may seem esoteric and pedantic. Its actually important, but its importance may not be clear to you until long after this class is over.
# 
# Dynamic typing means that you don't have to declare the type of a variable when you define it; python just figures it out based on how you are setting the variable. Lets say you set a variable. Sometime later you can just change the type of data assigned to a variable and python is perfectly happy about that. Since it won't be obvious until (possibly much) later why that's important, I'll let you marinate on that idea for a second. 
# 
# Here's an example of dynamic typing. What's the type of data assigned to voltage?

# <codecell>

type(voltage)

# <markdowncell>

# Lets assign a value of 2.7 (which is clearly a float) to voltage. What happens to the type?

# <codecell>

voltage = 2.7

# <codecell>

type(voltage)

# <markdowncell>

# You can even now assign a string to the variable voltage and python would be happy to comply.

# <codecell>

voltage = "2.7 volts"

# <codecell>

type(voltage)

# <markdowncell>

# I'll let you ruminate on the pros and cons of this construction while I change the value of voltage back to an int:

# <codecell>

voltage = 2

# <markdowncell>

# ## Coersion
# It is possible to coerce (a fancy and slightly menacing way to say "convert") certain types of data to other types. For example, its pretty straightforward to coerce numerical data to strings.

# <codecell>

voltageString = str(voltage)

# <codecell>

currentString = str(current)

# <codecell>

voltageString

# <codecell>

type(voltageString)

# <markdowncell>

# As you might imagine, you can go the other way in certain cases. Lets say you had numerical data in a string.

# <codecell>

resistanceString = "4.0"

# <codecell>

resistance = float(resistanceString)

# <codecell>

resistance

# <codecell>

type(resistance)

# <markdowncell>

# What would happen if you tried to coerce resistanceString to an int? What about coercing resistance to an int? Consider the following:

# <codecell>

resistanceString = "4.0 ohms"

# <markdowncell>

# Do you think you can coerce that string to a numerical type?
# ## On Being Precise with floats and ints
# Again, the following may seem esoteric and pedantic, but it is very important. So bear with me.
# Let's say you had some voltage data that looks like the following
# ``
# 0
# 0.5
# 1
# 1.5
# 2
# ``
# 
# Obviously, if you just assigned this data individually to a variable, you'd end up with the following types
# ``
# 0   -> int
# 0.5 -> float
# 1   -> int
# 1.5 -> float
# 2   -> int
# ``
# 
# But what if you wanted all of that data to be floats on its way in? You could assign the variable and then coerce it to type float:

# <codecell>

voltage = float(1)

# <markdowncell>

# But that's ugly. If you want what is otherwise an integer to be a float, just add a period at the end

# <codecell>

voltage = 1.

# <codecell>

type(voltage)

# <markdowncell>

# This point becomes important when we start operating on data in the next section.
# 
# ## Data Operations
# 
# What's the point of data if we aren't going to do something with it?  Let's get computing.

# <codecell>

a = 1

# <codecell>

b = 2

# <codecell>

c = a+b

# <codecell>

c

# <codecell>

type(a), type(b), type(c)

# <markdowncell>

# So we got a value of three for the sum, which also happens to be an integer. Any operation between two integers is another integer. Makes sense.
# 
# So what about the case where a is an integer and b is a float?

# <codecell>

a = 1

# <codecell>

b = 2.

# <codecell>

c = a + b

# <codecell>

c

# <codecell>

type(a), type(b), type(c)

# <markdowncell>

# You can do multiplication on numbers as well.

# <codecell>

a = 2

# <codecell>

b = 3

# <codecell>

c = a * b

# <codecell>

c

# <codecell>

type(a), type(b), type(c)

# <markdowncell>

# Also division.

# <codecell>

a = 1

# <codecell>

b = 2

# <codecell>

c = a / b

# <codecell>

c

# <markdowncell>

# **ZING!**
# 
# This is why type is important. Divding two integers returnes an integer: this operation calculates the quotient and floors the result to get the answer.
# 
# If everything was a float, the division is what you would expect.

# <codecell>

a = 1.

# <codecell>

b = 2.

# <codecell>

c = a / b

# <codecell>

c

# <codecell>

type(a), type(b), type(c)

# <markdowncell>

# There are operations that can be done with strings.

# <codecell>

firstName = "Johann"

# <codecell>

lastName = "Gambolputty"

# <markdowncell>

# When concatenating strings, we must explicitly use the concatenation operator +.  Computers don't understand context.

# <codecell>

fullName = firstName + lastName

# <codecell>

print fullName

# <codecell>

fullName = firstName + " " + lastName

# <codecell>

print fullName

# <markdowncell>

# There are other operations deined on string data. Use the *dir* comnand to find them. One example I'll show is the upper method. Lets take a look at the documentation.

# <codecell>

str.upper?

# <markdowncell>

# So we can use it to upper-caseify a string. 

# <codecell>

fullName.upper()

# <markdowncell>

# You have to use the parenthesis at the end because upper is a method of the string class.
# 
# For what its worth, you don't need to have a variable to use the upper() method, you could use it on the string itself.

# <codecell>

"Johann Gambolputty".upper()

# <markdowncell>

# What do you think should happen when you take upper of an int?  What about a string representation of an int?
# 
# That wraps up this lesson. We tried out the iPython shell and got some experience with ints, floats, and strings. Along the way we talked about some philosophy and how programming is like hammering.  
# 
# ## Miscellaneous scraps
# ## Pasting
# 
# You can paste things into the ipython console by copying text from your machine with **ctrl+c** and typing **%paste** at the iPython prompt.  The **%paste** is necessary syntax for multi-line clipboard deposits.
# 
# <codecell>
%paste
   
# <markdowncell>

