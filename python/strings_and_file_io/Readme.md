Python 5: Strings and File I/O
==============================

**Presented By : Tommy Guy**

Lesson goals:

1.  Examine the string class in greater detail.
2.  Use open() to open, read, and write to files.

Strings
-------

To start understanding the String class, let's use the built in help
system.

```python
  > help(str)
 Help on class str in module __builtin__:

 class str(basestring)
  |  str(object) -> string
  |  
  |  Return a nice string representation of the object.
  |  If the argument is a string, the return value is the same object.
  |  
  |  Method resolution order:
  |      str
  |      basestring
  |      object
  |  
  |  Methods defined here:
  |  
  |  __add__(...)
  |      x.__add__(y) <==> x+y
  |  
  |  __contains__(...)
  |      x.__contains__(y) <==> y in x
  |  
  |  __eq__(...)
  |      x.__eq__(y) <==> x==y
  |  
  |  __format__(...)
  |      S.__format__(format_spec) -> unicode
  |  
  |   __ge__(...)
  |      x.__ge__(y) <==> x>=y
  |  
  |  __getattribute__(...)
  |      x.__getattribute__('name') <==> x.name
  |  ]]
  |  __getitem__(...)
  |      x.__getitem__(y) <==> x[y]
  |  
  |  __getnewargs__(...)
  |  
  |  __getslice__(...)
  |      x.__getslice__(i, j) <==> x[i:j]
  |      
  |      Use of negative indices is not supported.
  |  
  ...
```

The help page for string is very long, and it may be easier to keep it open
in a browser window by going to the [online Python
documentation](http://docs.python.org/library/stdtypes.html#sequence-types-str-unicode-list-tuple-bytearray-buffer-xrange)
while we talk about its properties.

At its heart, a string is just a sequence of characters. Basic strings are
defined using single or double quotes.

```python
    > s = "This is a string."
    > s2 = 'This is another string that uses single quotes'
```

The reason for having two types of quotes to define a string is
emphasized in these examples:

```python
    > s = "Bob's mom called to say hello."
    > s = 'Bob's mom called to say hello.'
```

The second one should be an error: Python interprets it as `s = 'Bob'` then the
rest of the line breaks the language standard.

Characters in literal strings must come from the ASCII character set,
which is a set of 127 character codes that is used by all modern
programming languages and computers. Unfortunately, ASCII does not have
room for non-Roman characters like accents or Eastern scripts. Unicode
strings in Python are specified with a leading u:

```python
    > u = u'abcdé'
```

For the rest of this lecture, we will deal with ASCII strings, because
most scientific data that is stored as text is stored with ASCII.

Working with Strings
--------------------

Strings are iterables, which means many of the ideas from lists can also
be applied directly to string manipulation. For instance, characters can
be accessed individually or in sequences:

```python
    > s = 'abcdefghijklmnopqrstuvwxyz'
    > s[0]
    'a'
    > s[-1]
    'z'
    > s[1:4]
    'bcd'
```

They can also be compared using sort and equals.

```python
    > 'str1' == 'str2'
    False
    > 'str1' == 'str1'
    True
    > 'str1' < 'str2'
    True
```

In the help screen, which we looked at above, there are lots of
functions that look like this:

```python
    |  __add__(...)
    |      x.__add__(y) <==> x+y

    |  __le__(...)
    |      x.__le__(y) <==> x<y
```

These are special Python functions that interpret operations like \< and \+.
We'll talk more about these in the next lecture on Classes.

Some special functions introduce handy text functions.

**Hands on example**

Try each of the following functions on a few strings. What does the
function do?

```python
    > s = "This is a string"
    > s.startswith("This")
    > s.split(" ")
    > s.strip() # This won't change every string!
    > s.capitalize()
    > s.capwords()
    > s.lower()
    > s.upper()
```

File I/O
--------

Python has a built-in function called "open()" that can be used to
manipulate files. The help information for open is below:

    > help(open)
     Help on built-in function open in module __builtin__:

     open(...)
       open(name[, mode[, buffering]]) -> file object

       Open a file using the file() type, returns a file object.  This is the
       preferred way to open a file.

The main two parameters we'll need to worry about are the name of the
file and the mode, which determines whether we can read from or write to
the file. open returns a file object, acts like a pointer into the file.
An example will make this clear. In the code below, I've opened a file
that contains one line:

    (unix shell) $ cat testFile.txt
    abcde
    fghij

Now let's open this file in Python:

```python
    > fileHandle = open('testFile.txt','r')
```

The second input, 'r' means I want to open the file for reading only. I
can not write to this handle. The read() command will read a specified
number of bytes:

```python
    > s = fileHandle.read(3)
    > print s
    abc
```

We read the first three characters, where each character is a byte long.
We can see that the file handle points to the 4th byte (index number 3)
in the file:

```python
    > fileHandle.tell()
    3L
    > fileHandle.read(1)
    'd'
```

The file we are using is a long series of characters, but two of the
characters are new line characters. If we looked at the file in
sequence, it would look like "abcdenfghijn". Separating a file into
lines is popular enough that there are two ways to read whole lines in a
file. The first is to use the readlines() method:

```python
    > fileHandle.close() # close the old handle
    > fileHandle = open('testFile.txt','r')
    > lineArr = fileHandle.readlines()
    > lineArr
    ['abcde\n', 'fghij\n']
```

A very important point about the readline method is that it *keeps* the
newline character at the end of each line. You can use the strip()
method to get rid of the string.

File handles are also iterable, which means we can use them in for loops
or list extensions:

```python
    > f = open('testFile.txt','r')
    > l = [s.strip() for s in f]
    > l
    ['abcde', 'fghij']
    > f.close()
    > l = []
    > f = open('testFile.txt','r')
    > for s in f:
         l.append(s.strip())
```

These are equivalent operations. It's often best to handle a file one
line at a time, particularly when the file is so large it might not fit
in memory.

The other half of the story is writing output to files. We'll talk about
two techniques: writing to the shell and writing to files directly.

If your program only creates one stream of output, it's often a good
idea to write to the shell using the print function. There are several
advantages to this strategy, including the fact that it allows the user
to select where they want to store the output without worrying about any
command line flags. You can use "\>" to direct the output of your
program to a file or use "|" to pipe it to another program.

Sometimes, you need to direct your output directly to a file handle. For
instance, if your program produces two output streams, you may want to
assign two open file handles. Opening a file for reading simply requires
changing the second option from 'r' to 'w' or 'a'.

*Caution!* Opening a file with the 'w' option means start writing *at
the beginning*, which may overwrite old material. If you want to append
to the file without losing what is already there, open it with 'a'.

Writing to a file uses the write() command, which accepts a string.

```python
    > outFile = open('outfile.txt','w')
    > outFile.write('This is the first line!')
```

Another way to write to a file is to use writelines(), which accepts a
list of strings and writes them in order. *Caution!* writelines does not
append newlines. If you really want to write a newline at the end of
each string in the list, add it yourself.

Aside: The first exercise
=========================

Yesterday, we asked you to edit a file in place. Many of you asked how
this was possible. The answer is that it is not. You can use f.seek()
and f.tell() to verify that even if your file handle is pointing to the
middle of a file, write commands go to the end of the file in append
mode. The best way to change a file is to open a temporary file in
/tmp/, fill it, and then move it to overwrite the original. On large
clusters, /tmp/ is often local to each node, which means it reduces I/O
bottlenecks associated with writing large amounts of data.
