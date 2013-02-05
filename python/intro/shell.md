# Python 1 : The Shell, Variables, and Basic Data Types 

**Presented By: Katy Huff**

**Based on Lecture Materials By: Milad Fatenejad** 
**With help from Tommy Guy and Many More**

What is Python ?
================

Python is an interpreted (pre-compiled) language. Its simple, high
level, human readable interface speeds the programming process in
exchange for some computation time overhead. Python is implemented
mostly in the C programming language, so, as python develops, it is
increasingly possible to do everything in Python that is possible in C.
Python is also free and open source, so if you find a bug or generate a
useful module, the Python Software Foundation will likely be happy to
merge your changes into the language.

During this session you are going to learn about some very basics about
how to execute python code as well as some examples of the built-in
Python data types.

Built-in data types are the basic building blocks of Python programs.
They are really basic things like strings and numbers (either integers,
complex or floating point numbers). There are simple containers like
lists (think of lists as arrays or vectors), tuples and dictionaries.

Hello World
===========

First, we will use python ''interactively''. This means that we will
type commands directly into iPython. Once we start performing more
complicated tasks we will start writing Python scripts and programs in a
text editor, outside of the interpreter.

To get to the python shell, type **python** into the terminal.

```python
   >>> print "Hello World"
   Hello World
   >>> exit()
```

To get to the interactive python interpreter, a more sophisticated
python shell, type **ipython** into the terminal.

```python
   In [1]: print "Hello World"
   Hello World
   In [2]: exit
```

You can also put the commands in a **.py** file and execute that file in
the terminal by typing **python [filename]**

    $ gedit myfile.py &
    <edit myfile with the hello world program.>
    $ python myfile.py
    Hello World!

Pasting into iPython
====================

**Note:**

To paste text from another application (i.e. the internet browser) into
iPython:

1.  select text from the wiki
2.  copy with **ctrl+c**
3.  in iPython, type **%paste**

The code should paste and execute in iPython.

Variables
=========

Variables are names, while values are the data assigned to those names.

Questions : Variables and Values
--------------------------------

In the code snippet:

```python
    a=2
    b="string"
    c=a
```

 - What is the value of the variable `c`?
 - What is the value of the variable b ?
 - What is the name given to the variable 2 ?

(The last one is a trick, the value 2 has two names.)

Strings and Numbers
===================

It is really easy to make variables in python. For example, to create a
string, `s`, and print its value, simply type the following into
iPython:

```python
    s = "Hello World"
    print s
```

If you want to see what the type of a variable is, you can use the
built-in python function, `type`. Just enter

```python
   print type(s)
```

into iPython and you should see something like this:

```python
      <type 'str'>
```

This tells us that `s` is of type **str** (i.e. that `s` is a
string). Making numeric variables is equally easy and intuitive. Try
entering the following into IPython. Notice that the \# symbol is used
to start comments so everything after the pound sign is ignored.

```python
   i,r,c = -10, 3.5, 1.0 + 2j  # set i to -10, r to 3.5 and c to 1.0+2j
```

This one line sets the variable `i` to the integer -10 , `r` to the
floating point value 3.5 (a floating point number is just a
real/non-integer number) and `c` to the value 1.0 + 2j (Notice, how
easy and intuitive it is in python to set multiple variables to
something. You'll discover a lot of similar syntax that is designed to
make your life easier). Lets use the built-in type function to determine
the type of each of the three variables we just created:

```python
   print type(i), type(r), type(c)
```

This will give :
```python
    <type 'int'> <type 'float'> <type 'complex'>
```

This tells us that "i" is an integer, "r" is a floating point number,
and "c" is a complex number. As you can see, Python has built-in support
for imaginary numbers!

**Aside: Long integers** Another way python makes our lives easier is by
allowing integers to be arbitrary large. In languages like C/C++ and
FORTRAN integer variables can only store values up to a certain size.
But entering and manipulating the following forty digit number with
iPython is no problem:

```python
   i = 1234567890123456789012345678901234567890
   print i * 6
```

Operations in Python are defined by their type. For instance, look the
difference between these operations:

```python
   In[1]:  1 + 3
     4
   In[2]:  1.0 + 3
     4.0  # This is a float
   In[3]: "Hello " + "world"
     'Hello world'
   In[4]: 1 + "Hello"
   Traceback (most recent call last):
     File "<stdin>", line 1, in <module>
   TypeError: unsupported operand type(s) for +: 'int' and 'str'
```

In the first two cases, addition between numbers meant that 1 was added
to 3 using the standard type rules (float plus int = float). In the
third case, the command was string addition, which concatenates two
strings. The final case broke because an 'int' type can not be added to
a 'str' type. This is because it's unclear how to interpret an int as a
string: should it be the string representation, the ASCII character
code, or something else entirely?

One way to handle this is to explicitly convert the int into a string:

```python
     str(1) + "Hello"
```

Equivalent functions exist for converting to **int**, **float**, and
other types.

Basic data types in Python have a lot of functionality already built in.
For example, lets say that you are reading names from a file one line at
a time and that sometimes the names have leading and trailing spaces
that we want to strip away. We can just use the `strip` string method
to accomplish this. For example, type the following into iPython:

```python

  In[1]: name = "   Milad    "
  In[2]: print name + "is here"
        Milad     is here
```

Now enter `name.strip()` instead of `name`:

```python
  In[1]: print name.strip() + " is here"
   Milad is here
```

Notice that the extra spaces are gone. We used the `strip()` method,
which removes leading and trailing white space from strings. You can
think of a method as being a function that is attached to a particular
variable. You call methods by typing: `<variable>.<method name>`.

**Aside : Tab Completion**

Maybe you've noticed this already, but check out what happens you begin
typing a variable name (the first two letters of name, for example) and
press tab.

Convenient, right? This is also true of many built in functions.

Dynamic Typing
==============

Importantly, python is a **dynamically typed** language. That is, an
explicit type is not needed when creating a variable. Also, this means
that variables in Python which are initialized to a variable of one type
can be re-assigned to a variable of a different type. Try this:

```python
     sillystring = "What is the airspeed velocity of an unladen swallow?"
     print type(sillystring)
```

You'll see:

```python
     <type 'str'>
```

If you reassign silly string to an integer, what happens? That is, when
you type :

```python
    sillystring = 98    
    print type(sillystring)
```

You should see:

```python
     <type 'int'>
```

This is an interesting feature. Can you think of ways it can be helpful?
Are there ways it might be troublesome?

What is the type of sillystring be after this :

```python
    sillystring += 0.1
```

**Aside: In Place Equivalency**

What is the += syntax about? This is an in-place way to write `sillystring =
sillystring + 0.1`. It is common in a number of languages.

Importantly, though we do not explcity state them, variables always have
exactly one type. The number 98 is an **int**. For the variable holding
this value to be treated as a float, it must be assigned as **98.0**.

Questions : Dynamic Typing
--------------------------

Imagine that I first assign :

```python
    a=2
```

Then, I assign :

```python
    a="Welcome to the ministry of silly walks."
```

What has happened to the memory that was pointing to the number 2??

Getting Help
============

One of the really nice features in Python is that a lot of the help and
documentation is built into the code. Practically, this means that much
of the time you don't have to go digging through some web site to find
help. You can get help in Python using the `help` function. Lets look
at an example - enter

```python
    help(str.strip)
```

into IPython. You should then see documentation for the strip method pop
up. (NOTE: if you don't automatically return to the python interpreter,
just hit "`q`" to exit the help screen). You can also use the question
mark, "`?`", character to display the documentation as well. For
example, enter

```python
    str.strip?
```

into IPython to view the documentation.

Now try entering

```python
    help(str)
```

You should see documentation for the entire string type, including all
of the string methods. This can be useful when you are trying to perform
a specific task, but you don't know the right function to call. For
example, lets say we want to convert the string "cooper" to uppercase,
and we want to know if there is a string method which can do the job for
us. Start by typing "`help(str)`" to pull up the string documentation.
You can scroll through the string methods until you find a method called
"upper" which has documentation that looks like:

    |  upper(...)
    |      S.upper() -> string
    |      |      Return a copy of the string S converted to uppercase.

These lines tell us that the string class has a method called "upper"
which can be used to convert strings to uppercase. Now enter:

```python
    name = "cooper"   print name.upper()
```

At which point, you should see the word "COOPER" printed to the screen.

**Aside: Using Methods Directly on Data**

* * * * *

In the previous example, we first created a string variable, `name`,
assigned it the value "cooper", then used the `upper` string method to
obtain the uppercased version of the string. We didn't have to create a
variable, however. We could simply enter:

```python
    print "cooper".upper()
```

To generate the uppercased version.

As we saw above, the **str** type has a lot of documentation associated
with it, and we had to sift through most of it to find the upper method.
If we had a way to simply print all of the **str** methods, we could
have probably figured out that the `upper` method is what we wanted by
the name and in a lot less time. Luckily, python has a built in
function, "`dir`", for just this situation. The `dir` function takes
a type name and prints all of the methods associated. Try entering
"`print dir(str)`" to see a list of every method and variable
associated with the string class. You can ignore the methods that start
and end with double underscores for now. Try printing the methods
associated with the **int**, and **complex** types.

Finally, there are some really basic functions that are built right into
python that we have been using. For example, we used the "float" function
above to convert a string to a floating point number. You can see a list of
built in functions by entering `dir(__builtins__)`.  If you see something
interesting, such as the `zip` function, you can examine what it does using
help(zip).

Example : Manipulating Basic Data Types
---------------------------------------

Use the basic data types we've learned about along with the `help` and
`dir` functions to figure out how to do the following using either one
function or one method call:

- Take the absolute value of the number -1.4
- Begin with the string "a MaN and His DOG" and create the string "A man
  and his dog"
- Return the position of the character 'e' in the string "my test string"
  (The answer is 4, since `m` is is at position 0 not position 1)

