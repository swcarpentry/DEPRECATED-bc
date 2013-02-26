Python 9 : NumPy
================

**Presented By : Katy Huff**

**Based on Lecture Materials By Matthew Terry**

This class will cover a broad overview of NumPy with brief illustrative
examples geared towards getting people familiar with the basic use
cases. Since NumPy has many advanced features that may be useful to
experienced programmers, these notes will occasionally link to more
advanced examples that readers can peruse on their own time.

**Aside: Code Examples**

In all the examples below, we assume that import numpy has already been
executed. If any other modules are needed, we will import them
explicitly.

NumPy basics
============

What is NumPy ?
---------------

NumPy is a Python package implementing efficient collections of specific
types of data (generally numerical), similar to the standard array
module (but with many more features). NumPy arrays differ from lists and
tuples in that the data is contiguous in memory. A Python **list**, 
```[0, 1, 2]```, in contrast, is actually an array of pointers to Python
objects representing each number. This allows NumPy arrays to be
considerably faster for numerical operations than Python lists/tuples.

Creating NumPy Arrays
---------------------

Creating a NumPy array is as simple as passing a sequence to
numpy.array. You can also explicitly specify the data-type if the
automatically-chosen one would be unsuitable.

```python    

  >>> A = numpy.array([1, 2.3, 4])   
  >>> A.dtype 
  dtype('float64')   
  >>> B= numpy.array([1, 2.3, 4], dtype int)   
  >>> B.dtype   
  dtype('int32') 
```

As you might expect, creating a NumPy array this way can be slow, since
it must manually convert each element of a list into its equivalent C
type (int objects become C ints, etc). There are many other ways to
create NumPy arrays, such as **numpy.identity**, **numpy.zeros**,
**numpy.zeros\_like**, or by manually specifying the dimensions and type
of the array with the low-level creation function:

```python    

  numpy.ndarray((2, 3, 4), dtype=complex) # new 2x3x4 array of complex numbers 
```

For many of the examples below, we will be using **numpy.arange** which,
similar to the Python built-in function **range**, returns a NumPy array
of integers from 0 to N-1, inclusive. Like **range**, you can also
specify a starting value and a step.


```python

  >>> numpy.arange(2, 5)
  array([2, 3, 4])
  >>> numpy.arange(1, 5, 2)
  array([1, 3])
  >>> numpy.arange(1, 10, 2)
  array([1, 3, 5, 7, 9])
```

### Exercise : Create an Array

Create a NumPy array with values ranging from 0 to 10, in increments of
0.5. Don't forget to use **help()** to find useful functions!

Data types
----------

When creating a NumPy array, you supply a dtype ("data type"), or one is
chosen for you. There are a total of 21 different array scalar types,
which can be used to specify the dtype. In addition to the scalar type,
you may also specify byte order (little- or big-endian) or even multiple
scalar types to be used as a light-weight tuple (similar to a C struct).
For everyday use, however, you can just pass in the appropriate scalar
type and NumPy will figure it out. Some common scalar types include:

    int     # Python-compatible int (usually a C long)
    intc    # C int
    float   # Python-compatible float (usually a C double)
    single  # C float
    double  # C double
    complex # Python-compatible complex

List operations
---------------

For basic operations, NumPy arrays can be accessed just like Python
lists and tuples. This means that you can use the square brackets to
access elements, **len()** to access the size of the array, and so on.

```python

  >>> A = numpy.arange(5)
  >>> A
  array([0, 1, 2, 3, 4])
  >>> A[3]
  3
  >>> A[3] = 42
  >>> A
  array([ 0,  1,  2, 42,  4])
  >>> len(A)
  5
```

Arithmetic
----------

Since NumPy exists to perform efficient numerical operations in Python,
it stands to reason that NumPy arrays have all the usual arithmetic
operations available to them. These operations are performed
element-wise (i.e. the same operation is performed independently on each
element of the array).


```python

  >>> A = numpy.arange(5)
  >>> B = numpy.arange(5, 10)
  >>> A
  array([0, 1, 2, 3, 4])
  >>> B
  array([5, 6, 7, 8, 9])
  >>> A+B
  array([ 5,  7,  9, 11, 13])
  >>> B-A
  array([5, 5, 5, 5, 5])
  >>> A*B
  array([ 0,  6, 14, 24, 36])
```

In addition, if one of the arguments is a scalar, that value will be
applied to all the elements of the array.

```python
  
  >>> A = numpy.arange(5)
  >>> 2*A
  array([0, 2, 4, 6, 8])
  >>> A**2
  array([ 0,  1,  4,  9, 16])
  >>> A+10
  array([10, 11, 12, 13, 14])
```

Comparison
----------

Much like the basic arithmetic operations we discussed above, comparison
operations are perfomed element-wise. That is, rather than returning a
single boolean, comparison operators compare each element in both arrays
pairwise, and return an **array** of booleans (if the sizes of the input
arrays are incompatible, the comparison will simply return False). For
example:

```python

  >>> A = numpy.array([1, 2, 3, 4, 5])
  >>> B = numpy.array([1, 1, 3, 3, 5])
  >>> A == B
  array([ True, False,  True, False,  True], dtype=bool)
```

From here, you can use the methods .any() and .all() to return a single
boolean indicating whether any or all values in the array are True,
respectively.

Advanced Indexing
-----------------

In addition to the usual methods of indexing lists with an integer (or
with a series of colon-separated integers for a slice), NumPy allows you
to index arrays in a wide variety of different ways for more advanced
operations.

Multi-dimensional Indexing
--------------------------

Unlike Python lists and tuples, NumPy arrays can be multidimensional.
This complicates somewhat how they are indexed. To access a single
element, you simply pass in a comma-separated list of indices as a
subscript. However, there are many other things you can do when indexing
multidimensional arrays.

For instance, suppose you want the elements of a 2D array where the
first dimension is 1. NumPy makes this extremely simple:

```python

  >>> A = numpy.arange(16).reshape(4, 4)
  >>> A
  array([[ 0,  1,  2,  3],
         [ 4,  5,  6,  7],
         [ 8,  9, 10, 11],
         [12, 13, 14, 15]])
  >>> A[1]
  array([4, 5, 6, 7])
```

Now suppose you want the elements of that array where the **second**
dimension is 1. To do this, you can use "slices" of an entire dimension
as a placeholder by typing : as your first "index". Continuing from
above:

```python

  >>> A[:, 1]
  array([ 1,  5,  9, 13])
```

As you'd expect, this type of indexing can become quite complicated with
arrays of high dimension.

### Exercise : Selective Array Display

Using what we've learned about slicing and indexing, print just the
upper-left quarter of the array A above.

Indexing with Arrays
--------------------

NumPy arrays can be indexed with other arrays, using either an array of
indices, or an array of booleans of the same length. In the former case,
NumPy returns a view of the data in the specified indices as a new
array. In the latter, NumPy returns a view of the array with only the
elements where the index array is True. (We'll discuss the difference
between views and copies in a moment.) This makes normally-tedious
operations like clamping extremely simple.

Indexing with an array of indices:

```python

  >>> A = numpy.arange(5, 10)
  >>> A
  array([5, 6, 7, 8, 9])
  >>> A[[0, 2, 3]]
  array([5, 7, 8])
  >>> A[[0, 2, 3]] = 0
  >>> A
  array([0, 6, 0, 0, 9])
```

Indexing with an array of booleans:

```python
  
  >>> import random
  >>> A = numpy.array([random.randint(0, 10) for i in range(10)])
  >>> A
  array([10,  5,  1,  2,  3,  9,  3,  4,  9,  8])
  >>> A[A>5] = 5
  >>> A
  array([5, 5, 1, 2, 3, 5, 3, 4, 5, 5])
```

NumPy Gotchas
=============

NumPy has some important interface differences from Python lists and
tuples that can confuse new users of the library. Below are the most
notable of these.

Multiplication and Addition
---------------------------

As you may have noticed above, since NumPy arrays are modeled more
closely after vectors and matrices, multiplying by a scalar will
multiply each element of the array, whereas multiplying a list by a
scalar will repeat that list N times.

```python
  
  >>> numpy.arange(5)*2
  array([0, 2, 4, 6, 8])
  >>> range(5)*2
  [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
```

Similarly, when adding two NumPy arrays together, we get the vector sum
back, whereas when adding two lists together, we get the concatenation
back.

```python
  
  >>> numpy.arange(5) + numpy.arange(5)
  array([0, 2, 4, 6, 8])
  >>> range(5) + range(5)
  [0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
```

Views vs. Copies
----------------

In order to be as efficient as possible, NumPy uses "views" instead of
copies wherever possible. That is, NumPy arrays derived from another
base array generally refer to the ''exact same data'' as the base array.
The consequence of this is that modification of these derived arrays
will also modify the base array. You saw this above in how the result of
an array indexed by an array of indices is a ''copy'', but an array
indexed by an array of booleans is a ''view''. (Phew!)

Specifically, slices of arrays are always views, unlike slices of lists
or tuples, which are always copies.

```python

  >>> A = numpy.arange(5)
  >>> B = A[0:1]
  >>> B[0] = 42
  >>> A
  array([42,  1,  2,  3,  4])
  >>> >>> A = range(5)
  >>> B = A[0:1]
  >>> B[0] = 42
  >>> A
  [0, 1, 2, 3, 4]
```

### Exercise : Copy a NumPy Array

Figure out how to create a copy of a NumPy array. Remember: since NumPy
slices are views, you can't use the trick you'd use for Python lists,
i.e. copy = list[:].

Mathematics with NumPy
======================

Being designed for scientific computing, NumPy also contains a host of
common mathematical functions, including linear algebra functions, fast
Fourier transforms, and probability/statistics functions. While there
isn't space to go over ''all'' of these in detail, we will provide an
overview of the most common/essential of these.

Basics
------

All NumPy arrays have a collection of basic operations built-in. Most of
these can be used to operate only on a particular axis of the array, but
for simplicity, we will only show them in action on one-dimensional
arrays.

```python
  
  >>> import random
  >>> A = numpy.array([random.randint(0, 10) for i in range(10)])
  >>> A
  array([6, 9, 9, 4, 9, 8, 7, 9, 0, 3])
  >>> A.min()
  0
  >>> A.max()
  9
  >>> A.mean()
  6.4000000000000004
  >>> A.std() # standard deviation
  2.9732137494637012
  >>> A.sum()
  64
```

For 2-dimensional (or more) arrays, there are some other common
operations:

```python
  
  >>> A = numpy.arange(16).reshape(4, 4)
  >>> A
  array([[ 0,  1,  2,  3],
         [ 4,  5,  6,  7],
         [ 8,  9, 10, 11],
         [12, 13, 14, 15]])
  >>> A.T # transpose
  array([[ 0,  4,  8, 12],
         [ 1,  5,  9, 13],
       [ 2,  6, 10, 14],
       [ 3,  7, 11, 15]])
  >>> A.trace()
  30
```

There are many more methods like these available with NumPy arrays. Be
sure to consult the NumPy documentation before writing your own
versions!

Matrices
--------

So far, we've used two-dimensional arrays to represent matrix-like
objects. However, NumPy provides a specialized class for this. The
matrix class is almost identical to a two-dimensional NumPy array, but
has a few changes to the interface to simplify common linear algebraic
tasks. These are: \* The `*` operator is performs matrix multiplication
\* The `**` operator performs matrix exponentiation \* The property `.I`
(or the method `.getI()`) returns the matrix inverse \* The property
`.H` (or the method `.getH()`) returns the conjugate transpose

### Example: Solving a System of Linear Equations

```python

  >>> import numpy.linalg
  >>> A = numpy.matrix([[3, 2, -1], [2, -2, 4], [-1, .5, -1]])
  >>> B = numpy.array([1, -2, 0])
  >>> numpy.linalg.solve(A, B)
  array([ 1., -2., -2.])
```

Universal Functions
===================

Universal functions (also called ufuncs) are high-speed, element-wise
operations on NumPy arrays. They are, in essence, what allows you to
operate on NumPy arrays efficiently. There are a large number of
universal functions available covering most of the basic operations that
get performed on data, like addition, subtraction, logarithms, and so
on. Calling a ufunc is a simple matter:

```python

  >>> A = numpy.arange(1,10)
  >>> numpy.log10(A)
  array([ 0.        ,  0.30103   ,  0.47712125,  0.60205999,  0.69897   ,
          0.77815125,  0.84509804,  0.90308999,  0.95424251])
```

In addition to basic operation like above, ufuncs that take two input
arrays and return an output array can be used in more advanced ways.

ufuncs
------

### Exercise : Elementwise Operations

Using ufuncs, calculate the reciprocals of each element in the following
array:

```python
  
  [8.1, 1.6, 0.9, 4.3, 7.0, 7.3, 4.7, 8.2, 7.2, 3.0,
  1.4, 9.8, 5.7, 0.7, 8.7, 4.6, 8.8, 0.9, 4.4, 4.4]
```

External Resources
==================

NumPy has too many features to discuss here. However, there are plenty
of resources on the web that describe NumPy in detail. Here is a
selection of them:

 * [NumPy user's guide](http://docs.scipy.org/doc/numpy/user)
 * [NumPy Reference](http://docs.scipy.org/doc/numpy/reference/)
 * [NumPy For Matlab Users](http://www.scipy.org/NumPy_for_Matlab_Users)
 * [NumPy CookBook](http://www.scipy.org/Cookbook) 
