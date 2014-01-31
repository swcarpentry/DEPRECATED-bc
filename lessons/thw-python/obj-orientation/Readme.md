Python 6: Classes and Objects
=============================

**Materials orignially by Tommy Guy**

Basic Classes
-------------

The file manipulation example in the last lecture hid a pretty amazing
idea. When we opened a file using open(), it returned a new type we
hadn't seen before. The type had methods that we could use to access the
file all at once or one character at a time. It had other methods to
move around in the file, to close it, and to update it. It also had
data: what file is open and where in the file is the next readable byte.
The magic is that you, the programmer, don't have to think about the
details of the file implementation. You just have to use the methods
available to access the file. This thought process is the basis of
objects and object oriented programming.

Object oriented (OO) programming revolves around the create and
manipulation of objects that have attributes and can do things. They can
be as simple as a coordinate with x and y values or as complicated as a
dynamic webpage framework. They are a useful way to organize programs.
C++ and Java are OO languages. Even fortran is adding OO constructs in
newer standards. Here is the code for making a very simple class that
sets an attribute. Start a new file, call it myclass.py and type this
in.

```python
    class MyClass(object):
      def setA(self, A):
        self.A = A
```

Now, in the Python shell, lets import and use MyClass:

```python
    > import myclass
    > anObject = myclass.MyClass()  # The MyClass object is in the myclass module.
    > type(anObject)
    <class 'myclass.MyClass'>  # See, it's a new type!
    > anObject.A = 34          # Set the class variable A directly.
    > print anObject.A
    > anObject.setA('hello')   # Set the class variable A with the setter method.
    > print anObject.A
```

It will help to have a bit of object-oriented vocabulary to understand what just happened:
 - Class - user defined type (MyClass)
 - object - instance of a Class (tmp = MyClass(), tmp is an object)
 - method - a Class function, also called a member function (tmp.getA())
 - attribute - a Class variable (tmp.A)

Remember: you *write* a class and you *create* and object.

**Hands on Example**

Write an Atom class with functions set_number(number_string) and
area_code().

```python
    """Matrix defines a real, 2-d matrix."""

    class Matrix(object):
      """I am a matrix of real numbers"""

      def __init__(self,h,w):
          self._nrows = h
          self._ncols = w
          self._data = [0] * (self._nrows * self._ncols)

      def __str__(self):
          return "Matrix: " + str(self._nrows) + " by " + str(self._ncols)

      def setnrows(self, w):
           self._nrows = w
           self.reinit()

      def getnrows(self):
           return self._nrows

      def getncols(self):
           return self._ncols

      def reinit(self):
           self._data = [0] * (self._nrows * self._ncols)

      def setncols(self, h):
           self._ncols = h
           self.reinit()

      def setValue(self,i,j, value):
         if i < self._nrows and j < self._ncols:
             self._data[i * self._nrows + j] = value
         else:
             raise Exception("Out of range")

      def multiply(self, otherMatrix):
         ''' Perform matrix multiplication and return a new matrix.
         The new matrix is on the left. '''
         result = Matrix(self._nrows, otherMatrix.getncols())
         # Do multiplication...
         return result

      def inv(self):
         ''' Invert matrix '''
         if self._ncols != self._nrows: raise Exception("Only square matrices are invertible")
         invertedMatrix = Matrix(self._ncols, self._nrows)
         invertedMatrix.setncols(self._ncols)
         invertedMatrix.setnrows(self._ncols)
         # INVERT!
         return invertedMatrix
```

Note the "self" argument in all of the class methods. This is a pointer
to the current object. You have to use self to reference methods and
data in an object.

The Big Idea: Interface vs. Implementation
------------------------------------------

Users shouldn't have to know *how* your program works in order to use
it.

The interface is a *contract* saying what a class knows how to do. The
code above defines matrix multiplication, which means that
mat1.multiply(mat2) should always return the right answer. It turns out
there are many ways to multiply matrices, and there are whole Ph.Ds
written on performing efficient matrix inversion. The implementation is
the way in which the *contract* is carried out.

Constructors
------------

Usually you want to create an object with a set of initial values for
things. Perhaps an object needs certain information to be created. For
this you write a "constructor." In python, constructors are just methods
with a special name:

```python
    class MyClass(object):
      def __init__(self):
          ''' Initialize things '''
```

**Aside: Magic functions**

Methods with leading and trailing double underscores are "magic
functions" in python.

 - Iteration (for x in sequence) uses \_\_next\_\_
 - Slicing ie brackets) (a[1:2]) uses \_\_get\_\_
 - Calling ie parentheses (a(3)) uses \_\_call\_\_
 - Help uses \_\_doc\_\_

*Write an initializer for the Matrix class that sets the height and
width. Change the multiply and inv methods to use this compiler*

**Aside: Old vs New Style Classes**

It is worth noting that there are two types of classes in python: Old
style classes (OSC) and new style classes (NSC). NSC fix some conceptual
problems with OSC (typing, diamond inheritance, subclassing built in
types, etc). Consequently OSC are gone in python 3. This is not a cause
for concern or confusion as the difference are subtle and will not
affect you until you have written enough python to be comfortable with
the distinction. Or you can always subclass object and render the issue
moot. Below illustrates the fix in the typing system.

```python
    class OldStyleClass: # don't use this one
        def __init__(self):
            print "Old Style"

    class NewStyleClass(object): # <- use this one
        def __init__(self):
            print "New Style"

    ns = NewStyleClass()
    os = OldStyleClass()

    print type(os)
    print type(ns)
```

Class methods with variable numbers of arguments
------------------------------------------------

In the previous session you learned about the power of python functions.
The full power of functions (keyword arguments, variable length
arguments, etc) are available in classes. For converts from other
languages, be aware that python functions do not have a type signature
so function overloading is not available.

Subclassing
-----------

If you want a to create a Class that behaves mostly like another class,
you should not have to copy code. What you do is subclass and change the
things that need changing. When we created classes we were already
subclassing the built in python class "object."

For example, let's say you want to write a sparse matrix class, which
means that you don't explicitly store zero elements. You can create a
subclass of the Matrix class that redefines the matrix operations.

```python
    class SparseMatrix(Matrix):
    """I am a matrix of real numbers"""

    def __str__(self):
      return "SparseMatrix: " + str(self._nrows) + " by " + str(self._ncols)

    def reinit(self):
      self._data = {}

    def setValue(self,i,j, value):
       self._data[(i,j)] = value


    def multiply(self, otherMatrix):
       ''' Perform matrix multiplication and return a new matrix.
       The new matrix is on the left. '''
       result = SparseMatrix(self._nrows, otherMatrix.getncols())
       # Do multiplication...
       return result

    def inv(self):
       ''' Invert matrix '''
       if self._nrows != self._rcols: raise Exception("Only square matrices are invertible")
       invertedMatrix = SparseMatrix(self._ncols, self._nrows)
```

The SparseMatrix object is a Matrix but some methods are defined in the
*super class* Matrix.

Python guidelines for code formatting and pythonic conventions on class
behavior, naming, etc.  [python
conventions](http://www.python.org/dev/peps/pep-0008)
