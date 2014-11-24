---
layout: lesson
root: ../..
title: Python Reference
---

#### Basic Operations

*   Use `variable = value` to assign a value to a variable.
*   Use `print first, second, third` to display values.
*   Python counts from 0, not from 1.
*   `#` starts a comment.
*   Statements in a block must be indented (usually by four spaces).
*   `help(thing)` displays help.
*   `len(thing)` produces the length of a collection.
*   `[value1, value2, value3, ...]` creates a list.
*   `list_name[i]` selects the i'th value from a list.

#### Control Flow

*   Create a `for` loop to process elements in a collection one at a time:

        for variable in collection:
            ...body...

*   Create a conditional using `if`, `elif`, and `else`:

        if condition_1:
            ...body...
        elif condition_2:
            ...body...
        else:
            ...body...

*   Use `==` to test for equality.
*   `X and Y` is only true if both X and Y are true.
*   `X or Y` is true if either X or Y, or both, are true.
*   Use `assert condition, message` to check that something is true when the program is running.

#### Functions

*   `def name(...params...)` defines a new function.
*   `def name(param=default)` specifies a default value for a parameter.
*   Call a function using `name(...values...)`.

#### Libraries

*   Import a library into a program using `import libraryname`.
*   The `sys` library contains:
    *   `sys.argv`: the command-line arguments a program was run with.
    *   `sys.stdin`, `sys.stdout`: standard input and output.
*   `glob.glob(pattern)` returns a list of files whose names match a pattern.

#### Arrays

*   `import numpy` to load the NumPy library.
*   `array.shape` gives the shape of an array.
*   `array[x, y]` selects a single element from an array.
*   `low:high` specifies a slice including elements from `low` to `high-1`.
*   `array.mean()`, `array.max()`, and `array.min()` calculate simple statistics.
*   `array.mean(axis=0)` calculates statistics across the specified axis.

#### Community standards

As you begin to develop python packages (i.e., bundled collections of python code) that others are using, or that you're hoping other developers will contribute to, it's useful to adhere to python community standards. Some python community standards that you should be aware of (and ideally adhere to in your own python package) include:
*   [pep8](https://www.python.org/dev/peps/pep-0008): a style guide for python, which discusses topics such as how you should name variables, how you should use indentation in your code, and how you should structure your ``import`` statements, among many other things. Adhering to pep8 makes it easier for other python developers to read and understand your code, and to understand what their contributions should look like.
*   [numpydoc](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt): a standard for API documentation through docstrings used by numpy, scipy, and many other python scientific computing pacakges. Adhering to numpydoc helps ensure that users and developers will know how to use your python package, either for their own analyses or as a component of their own python packages. If you use numpydoc, you can also use existing tools to automatically generate HTML documentation for your API.
* [pypi](https://pypi.python.org/pypi) (the python package index) [pip](https://pypi.python.org/pypi/pip): standards for making your python package accessible and installable from the command line. Uploading releases of your python package to pypi and testing that they are installable with pip enables users to easily obtain working versions of your python packages, and developers to easily distribute their own tools that rely on your python package.
* [Semantic Versioning](http://semver.org/): a standard describing how to define versions of your python package. Using Semantic Versioning makes it easy for other developers to understand what is guaranteed to stay the same and what might change across versions of your python package. This in turn enables other developers to confidently build tools that depend on your python package.
