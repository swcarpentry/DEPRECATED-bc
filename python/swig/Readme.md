**Presented By Katy Huff**

**Based on Lecture Material by Kurt Smith**

[SWIG](http://www.swig.org) is a Simple Wrapper Interface Generator. SWIG parses
C/C++ header files and generates the appropriate wrapper code to make your C/C++
code directly callable from python (among other languages). This makes it a very
convenient way of exposing C/C++ code to python. For more detailed info, check
out SWIG's python
[documentation](http://www.swig.org/Doc1.3/Python.html#Python).

In C/C++, if you have the header files, you have everything you need to
link correctly. In theory, if you have the header files for some C/C++
code, you can generate the glue code to such that python can call that
code. All you need to do is parse the header, figure out the name
mangling, write the appropriate C code to make it python aware, and
write the appropriate python to load the C code correctly. To me this
sounds like a large undertaking. Thankfully, SWIG can do most of that
for you.

To use SWIG, you write a config file that tells SWIG which header files
to parse. When you invoke SWIG, you pass it the config file and the
desired output language. SWIG generates the appropriate C/C++ and python
glue code. You compile and link the generated C/C++ and import the
generated python module. Sounds simple enough, lets work through and
example.

**Note regarding C vs C++**

I realize that C and C++ are not the same language and this tutorial
will work with C++ and not C. C coders should console yourselves with
the a sense of smug superiority that there is no C analog to the
C++-python problems we will address. The techniques introduced should
still transfer toward making your C API more pythonic.

**Note regarding*nix vs Windows**

I have no experience with SWIG and Windows. I'm told it isn't
[difficult](http://www.swig.org/Doc1.3/Windows.html). All my programming and
building experience is in linux and I wasn't up to the challenge of figuring out
declspec just for this tutorial (extra credit for anyone that can
explain/justify declspec to me).

Illuminating Example
====================

Deciding what to have for breakfast is a computationally difficult
problem. Breakfast calculations taking excessively long have been known
to cause the phenomenon of "brunch." Since it is computationally
expensive, python is the the appropriate language for writing a
Breakfast library. Something closer to the metal, like C++, is. So we
write our library in C++ and then decide to generate python bindings
using SWIG.

To generate the python bindings to breakfast (called pyfast), we must
write a SWIG config file. This file names the python module that you
will be generating and tells SWIG which headers to parse and what
additional steps should be taken. For the sake of example, we have three
header files Bfast.hpp, Spam.hpp and Eggs.hpp, which define classes
Bfast, Spam and Eggs, respectively. Spam and Eggs inherit from Bfast.
Here is our initial pyfast.swg configuration file. (`` `*.i ``\` is also
a common extension for SWIG config files)

```c 
  %module pyfast
  %{
  #include "Bfast.hpp"
  #include "Eggs.hpp"
  #include "Spam.hpp"
  %}

  %include "std_string.i"

  %include "Bfast.hpp"
  %include "Eggs.hpp"
  %include "Spam.hpp"
```

The \#include's are normal preprocessor directives. The %include's are SWIG
directives that tell to actually generate the python bindings.  Notice that you
first \#include the file, then you tell swig to generate the bindings with the
%include directive. Lets start building things.  Below is a simple makefile for
%building our swig pyshapes module.  [Similar
%examples](http://www.dabeaz.com/cgi-bin/wiki.pl?SwigFaq/SharedLibraries) exist
%for other platforms.

```Makefile
    #!Lineno
    #!Makefile
    PYTHON_INCLUDE_DIR = "/usr/include/python2.6"

    pyfast:
            swig -python -c++ pyfast.swg
            g++ -fPIC -shared \
                    -o _pyfast.so \
                    food.cpp \
                    pyfast_wrap.cxx \
                    -I $(PYTHON_INCLUDE_DIR)
```

-   **Line 4**: We call swig on pyfast.swg, telling it that the headers
    are C++ (swig assumes C by default) and that we want to generate
    python bindings. This generates two files: pyfast.py (the python
    bindings) and pyfast_wrap.cxx (the C++ glue code).
-   **Line 5**: We build the shared library. -fPIC and -shared are
    compiler options for position independent code and building a shared
    library.
-   **Line 6**: Continuation of build command. We specify that the name
    of the shared library will be _pyfast.so. If you look in the
    pyfast.py file you will see that this is the assumed name for the
    shared library.
-   **Line 7**: Continuation of build command. We specify the source
    code being wrapped. If libbfast.so already exists, we can link to
    that rather than compiling in.
-   **Line 8**: Continuation of build command. We compile the actual
    glue code.
-   **Line 9**: Continuation of build command. We tell the compiler
    where to find the python header files. Depending on your platform
    configuration, you may need to explicitly tell the compiler where
    the python runtime library is.

If all went well, we should now have _pyfast.so and pyfast.py. Lets run
a test script to see how things worked. Download it here.

Python and C++ are Different Languages
======================================

Python and C++ are different languages (surprise!). They have different
conventions and different features. Don't be surprised if there isn't a
direct analog of some C++ feature in python (or vice versa).

Templates and SWIG
------------------

Metaprogramming is really nice in C++ because it allows you to write
general algorithms, but get specialized performance. It is nice that the
compiler generates templated code for you, but it is awkward when you
want to like to that generated code. SWIG only knows how to link python
to compiled object code. The way around this is that for every template
that occurs in the python bindings, you have to manually create an
instance. So if your API takes a vector of doubles, you have to make a
vector of doubles in the your swig config. Add the following code to
your swig config.

    %include "std_vector.i"
    %template std::vector<double>;

SWIG is aware of deque, list, map, pair, set, string, and vector. The [SWIG STL
documentation](http://www.swig.org/Doc1.3/Library.html#Library_stl_cpp_library)
is quite helpful.

Renaming
--------

You might have a function in C++ that shares a name with a python
keyword (print perhaps?). One solution around this problem is to rename
your function. You add the following code to your swig config:

    %rename('new_name') old_name

SWIG can be fairly aggressive when renaming things. The above code will
rename all functions named old_name (including class member functions).
You can make the renaming more specific by adding a function signature
and/or class resolution. Lets add the following renaming to our project.

    %rename("__str__") Bfast::string_rep() const;

Ignoring
--------

Sometimes you don't want to expose certain functionality to python.
Sometimes you can't get something to work and you want the problem to
Just Go Away. SWIG can ignore the problem

    %ignore lotsOfSpam(const Bfast&);

Lets put this into our project and see what happens.

You should be careful with %ignore's and %rename's as they tend to be
greedy. If you have class member functions that have the same signature,
they will get renamed/ignored. Also be careful with typedef's as SWIG
doesn't always know that (typedef int Int) Int's are int's.

## [Cross language polymorphism](http://www.swig.org/Doc1.3/Python.html#Python_directors) 

Typical SWIG wrapping consists of generating a proxy class in python
that handles dispatching calls to the compiled C++ library. This works
great so long as your interaction with the library is one way. Function
calls, class instantiations, normal things are all one way
communication. Using python to subclass C++ classes with virtual
functions requires two way communication. The feature you want to
investigate in this case is called "directors".

Auxillary code
--------------

SWIG has the capability of including code directly in the swig config
file.

[C/C++ Helper Functions](http://www.swig.org/Doc1.3/Python.html#Python_nn41)

    %inline

[Additional Python Code](http://www.swig.org/Doc1.3/Python.html#Python_nn42)

    %pythoncode

More Info
=========

The SWIG [C++ documentation](http://www.swig.org/Doc1.3/SWIGPlus.html) is quite
helpful. It discusses everything covered here in greater detail and has sections
specifically dealing with

-   [Swig and C++](http://www.swig.org/Doc1.3/SWIGPlus.html)
-   [Swig and Python](http://www.swig.org/Doc1.3/Python.html)
-   [Wrapping Overloaded Functions and
    Methods](http://www.swig.org/Doc1.3/SWIGPlus.html#SWIGPlus_overloaded_methods)
-   [Templates](http://www.swig.org/Doc1.3/SWIGPlus.html#SWIGPlus_nn30)
-   [Exception Handling](http://www.swig.org/Doc1.3/Customization.html#exception)
-   [STL Exceptions](http://www.swig.org/Doc1.3/Library.html#Library_stl_exceptions)

