---
title: Strings and File I/O in Python
---
# Edited by Jin

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
room for non-Roman characters like accents or Eastern scripts.

Unicode strings in Python are specified with a leading u:

```python
    > u = u'abcdé'
```

*Tip: Unicode and ASCII are different types of character sets that can be handled by programs.
*     For all intents and purposes, we will assume that strings are in ASCII (as per Python standard)
*     unless otherwise stated. Moreover, most scientific data is stored in ASCII format as well!


Working with Strings
--------------------

Strings are iterables, which means many of the ideas from lists can also
be applied directly to string manipulation. In essence, being 'iterable' allows strings
to be handled like the way lists are handled.


For instance, characters can be accessed individually or in sequences:

```python
    > s = 'abcdefghijklmnopqrstuvwxyz'
    > s[0]
    'a'
    > s[-1]
    'z'
    > s[1:4]
    'bcd'
```

They can also be compared using equals and sorted.

```python
    > 'str1' == 'str2'
    False
    > 'str1' == 'str1'
    True
    > sorted('grenade')
    ['a', 'd', 'e', 'e', 'g', 'n', 'r']
```

Since strings are iterables, you can loop through strings!
```python
    > s = 'A piece of string'
    > for characters in s:
    >>    print characters
          A
          
          p
          i
          e...
```

Like lists, we can also see if a character, or sequence of characters, is within our string.
In a list, we can test to see if a value exists: Python iterates through all elements of the list
and checks to see if our query is within the list. For instance,
```python
    > my_list = ['a', 1, 3]
    > 'a' in my_list
    True
    > 4 in my_list
    False
```

For strings, we rely on a very similar syntax to check if our character(s) are within the string:
```python
    > s = 'The cat in the hat'
    > 'T' in s
    True
    > 'cat' in s
    True
    > 'Hat' in s 
    False
```
*Caution*
Note that in the last case, although the word 'hat' is clearly within the string, the searches are case-sensitive!
When we use these methods, make sure to search through carefully.


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

Python has a built-in function called `open()` that can be used to
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
that contains two lines:

    (unix shell) $ cat testFile.txt
    abcde
    fghij

Now let's open this file in Python:

```python
    > fileHandle = open('testFile.txt','r')
```

The second input, 'r' means I want to open the file for reading only. I
can not write to this handle. The `read()` command, associated to this fileHandle object, 
will read a specified number of bytes and output, interestingly, a string!

```python
    > s = fileHandle.read(3)
    > print s
    abc
```
By calling `read()` from the fileHandle object, we read the first three characters, where each character is a byte long.
The `read()` function iterates through each character and churns out a string from the file to screen.

*Caution*: Note that using read() will move along the file stream and the next time we call read(), we'd point to the
           (n+1)th byte. For example, if we invoke read(3), we have read 3 bytes, and the next time we call read(), we
           will read from the 4th byte and onward.

In fact, we can see that the file handle points to the 4th byte (index number 3; 
remember that iterables start at index 0?)

```python
    > fileHandle.tell()
    3L
    > fileHandle.read(1)
    'd'
```

The file we are using is a long series of characters, but two of the
characters are new line characters. If we looked at the file in
sequence, it would look like "abcde\nfghij\n". Separating a file into
lines is popular enough that there are two ways to read whole lines in a
file. The first is to use the `readlines()` method:

```python
    > fileHandle.close() # close the old handle
    > fileHandle = open('testFile.txt','r')
    > lineArr = fileHandle.readlines()
    > lineArr
    ['abcde\n', 'fghij\n']
```
A very important point about the readlines method is that it *keeps* the
newline character, `\n`, at the end of each line. You can use the `strip()`
method to get rid of the newline character.

*Caution*: confusingly enough, Python has both a readline() and readlines() function. Readlines(), as shown above,
           will read the entire file and separate each line to an individual element of a list. On the other hand,
           readline() will read one line at a time, and like the command read(), it will keep a pointer to the next line.

Example:

    (unix shell) $ cat testFile.txt
    abcde
    fghij

```python
    > file = open('testFile.txt')
    > file.readlines()
    ['abcde\n', 'fghij\n']
    > file.close()
    
    > sameFile = open('testFile.txt')
    > sameFile.readline()
    abcde
    > sameFile.readline()
    fghij

```

Notice that when we called `readline()` the second time for sameFile, it outputs the second line of the file?

File handles are also iterable, which means we can use them in for loops or list extensions:

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
assign two open file handles. Opening a file for reading # writing?
simply requires
changing the second option from 'r' to 'w' or 'a'.

*Caution!* Opening a file with the 'w' option means start writing *at
the beginning*, which may overwrite old material. If you want to append
to the file without losing what is already there, open it with 'a'. Writing to
a file uses the `write()` command, which accepts a string. Suppose that we have our
friendly file, testFile.txt:

    (unix shell) $ cat testFile.txt
      abcde
      fghij

If we open this file using the 'w' mode, then it will overwrite from the beginning!

```python
    > outFile = open('testFile.txt','w')
    > outFile.write('This is the first line!\n')
    > outFile.close()
    
```

So if we read this file now, it'd be...

    (unix shell) $ cat testFile.txt
      This is the first line!

Now let's open the file in 'a' mode this time.
```python
    > outFile = open('testFile.txt','a')
    > outFile.write('This is the second line!\n')
    > outFile.close()
    
```
Because we wrote to this file in append mode, we get:

    (unix shell) $ cat testFile.txt
      This is the first line!
      This is the second line!

Aside: The first exercise
=========================

Yesterday, we asked you to edit a file in place. Many of you asked how
this was possible. In a unix shell, this is possible by using the command `sed`, but
this is beyond the scope of our lesson - also, this is not a Python command!

You can read more about `sed` [here](http://www.grymoire.com/Unix/Sed.html),
and I'll put a short piece of code on how to use sed at the bottom.

*For the super keen*




You can use `seek()` and `tell()` to verify that even if your file handle is pointing to the
middle of a file, write commands go to the end of the file in append mode. Sometimes, if space is
not an issue, I would suggest making a copy of the file on-the-fly in Python. However, if it's a big
file (e.g. >1000 lines), you can read it line-by-line, making the file I/O very efficient.

Example: 

First, count the number of lines in a unix shell by using `wc -l <FILENAME>`. `wc` is a word-count command in unix,
and the `-l` option specifies line counting.

    (unix shell) $ wc -l bigFile.txt
      10000 bigFile.txt
  
bigFile.txt has 10000 lines! But we want to make some changes on the fly and keep a copy of bigFile.txt. An easy
way to do it is to keep one file handle for bigFile.txt but to create another file handle for the new file. Alternatively,
to reduce the hold on memory, you can keep a temporary copy of the file in your disk. 

I'll show you the on-the-fly method where we add the string " NOW EDITED " at the end of each line.

```python
    > fileHandle = open('bigFile.txt', 'r')
    > newFile = open('bigFile_Copy.txt', 'w')
    > for lines in fileHandle:
        lines = lines + " NOW EDITED \n"
        newFile.write(lines)
    
    > fileHandle.close()
    > newFile.close()

```

Aside: The sed method
=======================

`sed` is a built-in unix command for manipulating text files.
If we wanted to change a file in-place, then we can use substitute:

    (unix shell) $ sed -i 's/hello/bye/g' file.txt
  
In this case, for every instance of the word 'hello' in file.txt, sed will replace it with 'bye'.
By using `-i`, we invoke an in-place editing option, with `s/` asking for substitution. The syntax is:

    (unix shell) $ sed -i 's/<Pattern to match>/<New pattern to use>/g' [FILENAME]

With sed, we can also append lines. Note, for sed appending, we must put these commands as separate lines:

    (unix shell) $ sed '/<Pattern to Match>/ a\
                  <New Pattern>
                  ' [FILENAME]

So at first, we give sed the pattern to match, and the `a\` command to flag appending. Then in a new line, we put our
new pattern of interest. Once we submit, we close with a closing quotation mark, and the filename to append to.

    (unix shell) $ cat newFile.txt
      abcde
      fghij
      klmnop
      
                   sed '/fghij/ a\
                   HELLO WORLD!!
                   ' newFile.txt
                   
                   
                   cat newFile.txt
      abcde
      fghij
      HELLO WORLD!!
      klmnop

        








