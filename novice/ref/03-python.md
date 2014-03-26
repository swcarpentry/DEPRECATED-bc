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
