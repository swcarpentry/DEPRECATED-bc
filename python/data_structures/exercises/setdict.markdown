# Sets and Dictionaries

Warmup by Greg Wilson. Exercises by Orion Buske, with thanks to Tommy Guy and Jon Pipitone.

---

## Warmup

The following dictionary maps lowercase letters in the English alphabet to Morse
code symbols:
    
```python
morse_code = {'a':'.-', 'b':'-...', 'c':'-.-.', 'd':'-..', 'e':'.', 'f':'..-.', 
'g':'--.', 'h':'....', 'i':'..', 'j':'.---', 'k':'-.-', 'l':'.-..', 'm':'--', 
'n':'-.', 'o':'---', 'p':'.--.', 'q':'--.-', 'r':'.-.', 's':'...', 't':'-',
'u':'..-', 'v':'...-', 'w':'.--', 'x':'-..-', 'y':'-.--', 'z':'--..'}
```

### 1. 
Write a function, `encode`, that encodes a lowercase English word in Morse code,
leaving a space between each encoded character, e.g.:

```python
>>> encode("sos")
'... --- ...'
```

**Answer:**
<span class="answer">
Here's the straightforward long version:
```python``
def encode(word): 
    # we'll store the encoded characters in a list
    encoded = []
    for c in word: 
        encoded.append(morse_code[c])
    
    # and then join the characters with spaces
    return " ".join(encoded)
```

Here's a shorter version that uses a fast, efficient, list comprehension that we
won't go into detail describing, but you can google around to find out more:

```python
def encode(word): 
    encoded = [morse_code[c] for c in word]
    return " ".join(encoded)
```
</span>


### 2. 
What happens to characters that don't appear in the dictionary? Modify the
function to leave those characters untouched.

**Hint:**
<span class="hint">
```python
help(dict.get)
```
</span>
 
**Answer:**
<span class="answer">
```python
def encode(word):
    encoded = []
    for c in word: 
        encoded.append(morse_code.get(c,c))

    return " ".join(encoded)
```

The short version can similarly be modified.
</span>


### 3. 
Write a function, `decode`, to decode a string of morse code, as produced by the
encode function, to the english alphabet.

**Hint:**
<span class="hint">
You'll need to 'invert' the `morse_code` dictionary so that you can map from
morse code alphabet to the english alphabet. 

**Hint:**
<span class="hint">
```python
    help(dict)
```
</span>
     
**Hint:**
<span class="hint">
```python
    help(str.split)
```
</span>

**Answer:**
<span class="answer">
String appending is very slow in python, so we append to a list and then join on
an empty string:

```python
# invert the dictionary
english = {}
for (k,v) in morse_code.items(): english[v] = k

def decode(morse_word):
    decoded = []
    characters = morse_word.split()

    for c in characters: 
        decoded.append(english[c])
    
    return ''.join(decoded)
```

As before, this can also be easily implemented with a list comprehension:

```python
def decode(morse_word):
    decoded = [english[c] for c in morse_word.split()]
    return ''.join(decoded)
```
</span>


### 4. 

In what situation would inverting the dictionary as we have done above lead to a
loss of data?

**Answer:**
<span class="answer">
When there are multiple keys that point to the same value.  Since values
become keys and keys become values unless care is taken to somehow combine the
original keys into a single new value then only one of the original keys will
remain as a value.  
</span>


## Exercises

Danielle is a developer and was hired by a new tech startup. The company makes
mobile devices, and they have been loading them with Google's Android operating
system. However, they hired Danielle and a few other developers to try to do
better. They've decided to create a new, light-weight Linux distribution
tailored to the needs of the devices the company makes. The project is going
well, but they desperately need a package manager to deal with all the packages
that the company might want to put on their devices. Needless to say, Danielle
has been put in charge of this, and we're going to help her.


Many operating systems distribute software as a unit called a _package_. Since
software requires other software in order to run, packages list other packages
as _dependencies_. Since we don't want the actual computer user to have to track
and install all these dependencies, we have package managers. We tell the
package manager we want to install a certain package, like `pytables` (a Python
package), and it finds (and installs) all the dependencies that aren't already
installed on the system (like `numpy`, `numexpr`, `cython`, etc.), and then
installs `pytables`. Since those dependencies have dependencies, which in turn
have their own dependencies, this can quickly get hairy.


Luckily, they don't need anything super fancy. In particular, we are going to
completely ignore package version numbers. They should also probably use a
simple database to store and access this information, but for now they're going
to use a simple line-based text file. This makes it easy to view, edit, and
version control, but isn't as scalable or manageable. You can find a sample
dependency file (aside: these are all real software dependencies
from a few recent projects), in [setdict_dependencies.txt](setdict_dependencies.txt)  


### Exercise 1: reading in the file


Each line of the dependency file records a single dependency:
    
    package_name[TAB]required_package_name

Comment lines begin with the pound character ('#') and should be ignored. They
chose this file format instead of having a single line per package that listed
_all_ the dependencies because they had used subversion before and knew that
multiple users changing the same line can cause headaches, but changing
different lines is usually fine.

1.  Please create a `setdict` folder in your personal folder in the course
repository to submit your responses. For this exercise, we'll be writing a file
called `package_manager.py` that should be put in this folder and committed to
the repository.

2.  We would like to be able to quickly report all the dependencies of a given
package, so a dictionary seems like a good way to store this information. Add a
function called `read_dependencies` in `package_manager.py` that accepts an
opened dependency file in the above format and returns a dictionary that maps a
package name to the set of packages it directly requires. The set should contain
just first-level dependencies, so `Y` should be in `dependencies[X]` if and only
if `X[TAB]Y` is line in the file. Thus, the following example should hold:

   ```python    
       >>> with open("dependencies.txt") as ifp:
       ...    dependencies = read_dependencies(ifp)
       ...
       >>> print sorted(dependencies["genomedata"])
       ['hdf5', 'numpy', 'pytables', 'python']
   ``` 
   
   **Hint:**
   <span class="hint">
   Some pseudocode:
   ```
   def read_dependencies(reader):
       initialize dictionary
       for line in reader:
           parse line into package and dependency
           add dependency to set of dependencies for package
       return dictionary
   ```
   </span>


3.  If you haven't seen the "with" keyword before, here is a simple description.
It was added in Python 2.6 and provides a nice way of simplifying the following
pattern:

   ```python
       open(file)
       ...
       close(file)
   ```
   as
   ```python
       with open(file) as variable:
          ...
   ```

   Not only does the second form clarify the block of code that uses the file,
   it will automatically close the file when Python exits the code block. Under
   what situation does the second form close opened files that the first form
   does not?

   
   **Answer:**
   <span class="answer">
   If an exception is raised in the `...` block, Python will bail out until it
   finds a `try/catch` to catch the exception. In the first form, the file
   will remain open, which we almost certainly won't want. The second form
   will still nicely close the file for us, even in this case.
   </span>

4.  The second thing in the above code that might seem strange is that we
printed the sorted output. Why did we do this instead of just printing the set
returned by `dependencies["genomedata"]`?

   **Answer:**
   <span class="answer">
   Because sets aren't ordered, there is no way to predict what would be
   displayed from printing the set directly. Printed the sorted elements of
   the set ensures one can compare the outputs easily and is necessary if the
   example were to be used in some sort of unit test.
   </span>


### Exercise 2: finding dependencies
The core functionality of a package manager is that given a package the user
wants to install, it returns a set of other packages that need to be installed.
Normally, these other packages need to be installed _first_, so the problem ends
up being harder, but we'll simplify the problem slightly and say that all the
dependencies just need to be installed before the package can run. Thus, we
don't need to worry about tracking the order dependencies need to be installed
in, just which ones need to be installed. (This idilic case is like assuming all
the dependencies use runtime libraries that get installed in a default
location.)


1.  Add another function, `get_dependencies`, to `package_manager.py` that takes
two arguments, a dictionary mapping packages to their first-level dependencies
and the name of a package, and returns a set of **all** the dependencies of that
package. Not just the first-level dependencies, but also _those dependencies'
dependencies_, etc.

   You can assume there are no cycles (X requires Y, Y requires Z, Z requires X)
   or co-dependencies (X requires Y, Y requires X).

   This is essentially tree traversal, so if you're unfamiliar with tree
   traversal, check out the [Wikipedia
   page](http://en.wikipedia.org/wiki/Tree_traversal). The problem naturally
   lends itself to a solution that is recursive (having a function that calls
   itself). Once you look at the given package's dependencies, you need to look
   at the dependencies of each of those packages. Wouldn't it be awesome if we
   had a function that, given a package, returned a set of its dependencies? If
   so, we could just call it on each dependency and union the results. Oh, wait!
   We do! It's the function we're writing...


   Here's some expected output:
    ```python 
    >>> required_packages = get_dependencies(dependencies, "segway")
    >>> assert "lsf drmaa" in required_packages
    >>> print sorted(required_packages)
    ['cython', 'drmaa', 'gcc', 'gfortran', 'gmtk', 'hdf5',
    'lapack', 'lsf drmaa', 'numexpr', 'numpy', 'pytables',
    'python', 'python-nose', 'sge', 'zlib']
    ```

   **Hint:**
   <span class="hint">
   Some pseudocode to help...
   ```
   def get_dependencies(dependencies, package):
       look up direct dependencies of package
      if there aren't any:
        return an empty set
      else:
        initialize a set with the direct dependencies of package
        for every direct dependency, x, of package:
          add to the set all the items in get_dependencies(dependencies, x)
        return the set
   ```
   </span>

2.  In this simple recursive implementation, you explore the full dependency
tree. Unfortunately, this isn't very efficient. For example, `segtools` requires
`genomedata` and `pytables`, `genomedata` also requires `pytables`, and
`pytables` requires a number of other packages. You only need to explore the
dependency tree of `pytables` once. However, we're exploring both branches
completely separately, so that work is repeated! But hold on, since it is easier
to get the simpler, less efficient case correct, we're going to use it to test
our revised, more efficient version.

   Rename `get_dependencies` to `get_dependencies_slowly` and create a copy
   called `get_dependencies_quickly`. Tweak `get_dependencies_quickly` to avoid
   this redundant work.


   **Hint:**
   <span class="hint">
   The problem with the original function was that `get_dependencies` was
   called for _every_ dependency, even if that dependency was already
   explored. We need to keep track of which dependencies have already been
   explored, and we shouldn't call `get_dependencies` on those that already
   have been.
   </span>

   **Hint**
   <span class="hint">
   This can be implemented by passing around more data, a set of dependencies
   that have already been explored. You can give a default value to this
   argument so that it is initialized correctly on the first call, e.g.:
   
   ```python
   def get_dependencies_quickly(dependencies, package, _explored=None):
       if _explored is None:
           _explored = set()
       ...
   ```
   
   We do this `None` trickery with _explored instead of having
   `_explored=set()` as a default value for a very subtle reason. Basically,
   Python only evaluates the default values of arguments once, so if we add
   items to _explored, then subsequent calls to `get_dependencies_quickly`
   will default to _explored being a set with things already in it, rather
   than an empty set as you might hope. Thanks to David Boom for pointing this
   out.
   
   You might not have seen variables that start with an underscore ("_")
   before. It is the Python convention to start the name of private, internal
   variables with an underscore. It doesn't actually change the visibility of
   the variable, it just lets other programmers know that they probably
   shouldn't mess with it.
   </span>


3. Add an `if __name__ == '__main__':` block to `package_manager.py` which loads
a dependency file given as the first command-line argument (`sys.argv[1]`),
loads in the dependencies, calls `get_dependencies_slowly` and
`get_dependencies_quickly` for every package found (every key in the
dictionary), and asserts that the two results are the same for every package.
This is by no means an ideal way to test this, but it will suffice for our
purposes. Run code>package_manager.py` on the dependency file we gave a link to
at the beginning of the exercises with something like:  

  ```bash
  $ python package_manager.py /path/to/dependencies.txt
  ```

4.  Since we don't really want to expose the users of our package to our two
recursive functions, we should have `get_dependencies` somehow refer to
`get_dependencies_quickly`. There are two options for how we can implement this,
one with assignment and one by adding a level of abstraction:
    
   ```python
    get_dependencies = get_dependencies_quickly
   ```

   or

   ```python 
    def get_dependencies(dependencies, package):
        return get_dependencies_quickly(dependencies, package)
   ```
    
   Which should we choose, and why?

   **Answer:**
   <span class="answer">
   You might be tempted to just use assignment. It is simpler and slightly
   more efficient, but then if a user imports `get_dependencies` and then
   types `help(get_dependencies)`, they will see the extra argument that we
   use for keeping track of explored packages while recursing. In this case,
   adding the extra level of abstraction provides us more freedom and more
   clearly separates our implementation from our interface, so this is what we
   did.
   </span>

5. One final wrinkle is that users aren't really interested in _every_
dependency of a new package. Rather, they tend to be interested in which
dependencies they need to install. Change `get_dependencies` to allow an
optional argument, `installed`, that is a set of some packages that are already
installed. To make it easier for the user, the `installed` set does not have to
contain all installed programs, and you should assume the dependencies for those
programs are satisfied.

   Thus, because having `pytables` installed means all the dependencies for
   `genomedata` are already installed:
   ```python
   >>> required_packages = get_dependencies(dependencies, "genomedata",
   ...                                      installed=set(["pytables"]))
   >>> print sorted(requried_packages)
   []
   ```
    
   You might be tempted to pass `installed` as the parameter you added to
   `get_dependencies_quickly`, unfortunately this is not quite correct. Why
   isn't it correct to just skip packages that are already installed?

   **Hint:**
   <span class="hint">
   What happens if package A requires B and C, both of which require D, and B is
   already installed? Should D be in the set of packages that need to be
   installed?
   
   There are two relatively straightforward ways of solving this problem, one
   of which involves passing a set into `get_dependencies_quickly`, and one
   that doesn't. Neither is wrong, one is more efficient, and one might be
   considered better style. What are the two options, which did you choose,
   and why? Write your response in a comment in `get_dependencies` (not in the
   docstring, since we don't want `help(get_dependencies)` to print this), and
   then see our opinion below.
   </span>

   **Our opinion:**
   <span class="answer">
   It is more efficient to first find the union of all the dependencies of the
   installed packages, and then pass this as the "already-explored" argument
   to `get_dependencies_quickly`. However, in our implementation, the
   "already-explored" argument is only there to allow efficient recursion
   (which is why we named it with an underscore), so we want to avoid taking
   advantage of a particular aspect of the implementation that might change
   later. It is a hard call, but in this case, we went with the less-efficient
   code that doesn't take advantage of the function's extra argument.
   </span>

6. Your program should now work, and you should probably have all the sample
   outputs in your '__main__' block with `assert` statements verifying your
   output is correct (or unit tests with `nose` if you're exceptional). To allow
   you to improve your Python coding style, below is our code for the
   `read_dependencies` function. Please look over it carefully and note
   situations in which this code is more (or less) clear than the code you
   originally wrote. If you think you have a better way of solving it, please
   post on the forum. Then, using what you learned, rewrite your
   `get_dependencies_quickly` function with better style.  Leave
   `get_dependencies_slowly` alone since 1) you know it is correct and this lets
   you test your revised code and 2) it lets us see how much you improved.

   Our `read_dependencies` code:
   ```python
    def read_dependencies(reader):
        deps = {}
        for line in reader:
            line = line.strip()
            if not line or line.startswith('#'): continue
            # Ensures exactly two fields
            package, dependency = line.split('\t')
            if package not in deps:
                deps[package] = set()
            deps[package].add(dependency)
    
        return deps
   ```

