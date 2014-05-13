---
layout: lesson
root: ../..
---




# Creating Functions

If we only had one data set to analyze, it would probably be faster to load the file into a spreadsheet and use that to plot some simple statistics. 
But we have twelve files to check, and may have more in future. In this lesson, we'll learn how to write a function so that we can repeat several operations with a single command.

## Objectives

* Define a function that takes arguments (parameters).
* Return a value from a function.
* Test and debug a function.
* Explain what a call stack is, and trace changes to the call stack as functions are called.
* Set default values for function parameters.
* Explain why we should divide programs into small, single-purpose functions.
* Defining a function

Let's start by defining a function `fahr_to_kelvin()` that converts temperatures from Fahrenheit to Kelvin:


```r
fahr_to_kelvin <- function(temp) {
    kelvin <- ((temp - 32) * (5/9)) + 273.15
    kelvin
}
```


The definition opens with the name of your new function, which is followed by the call to make it a `function()` and a parenthesized list of arguments. You can have as many input parameters as you would like (but too many might be bad style). The body, or implementation, is surrounded by curly braces `{ }`. In many languages, the body of the function - the statements that are executed when it runs - must be indented, typically using 4 spaces. While this is not a mandatory requirement in R coding, we strongly recommend you to adopt this as good practice.

When we call the function, the values we pass to it are assigned to those variables so that we can use them inside the function. The last line within the function is what R will evaluate as a returning value. Remember that the last line has to be a command that will print to the screen, and not an object definition, otherwise the function will return nothing - it will work, but will provide no output. For example, let's try running our function. Calling our own function is no different from calling any other function:


```r
fahr_to_kelvin(32)
```

```
## [1] 273.1
```

```r
paste("boiling point of water:", fahr_to_kelvin(212))
```

```
## [1] "boiling point of water: 373.15"
```


We've successfully called the function that we defined, and we have access to the value that we returned. However, it the function was redefined as follows


```r
fahr_to_kelvin <- function(temp) {
    kelvin <- ((temp - 32) * (5/9)) + 273.15
}
```


Now typing


```r
fahr_to_kelvin(32)
```


Will return nothing. Why?

--> In Python lessons, now would come debugging…
==============================


## Composing Functions

Now that we've seen how to turn Fahrenheit into Kelvin, it's easy to turn Kelvin into Celsius:


```r
kelvin_to_celsius <- function(temp) {
    celsius <- temp - 273.15
    celsius
}

paste("absolute zero in Celsius:", kelvin_to_celsius(0))
```

```
## [1] "absolute zero in Celsius: -273.15"
```


What about converting Fahrenheit to Celsius? We could write out the formula, but we don't need to. Instead, we can compose the two functions we have already created:


```r
fahr_to_celsius <- function(temp) {
    temp_k <- fahr_to_kelvin(temp)
    result <- kelvin_to_celsius(temp_k)
    result
}

paste("freezing point of water in Celsius:", fahr_to_celsius(32))
```

```
## [1] "freezing point of water in Celsius: 0"
```


This is our first taste of how larger programs are built: we define basic operations, then combine them in ever-large chunks to get the effect we want. 
Real-life functions will usually be larger than the ones shown here - typically half a dozen to a few dozen lines - but they shouldn't ever be much longer than that, or the next person who reads it won't be able to understand what's going on. __Modular programming__


### Challenges

As we've seen in our print statements, we can use `paste` to concatenate strings, `paste(a, b, sep = "")` is `ab`. __Note__: the `sep` can be an important value to define! What is the default? What can `sep` be?

1. Write a function called `fence` that takes two arguments called `original` and `wrapper` and returns a new string that has the `wrapper` character at the beginning and end of the `original`:





```r
fence("name", "*")
```

```
## [1] "*name*"
```


## String splits

If the variable s refers to a string, then we can parse the string into its separate components - each of the characters. Base R has a function called `strsplit` that can be used to break up strings, into smaller chunks. 


```r
pangram <- "the quick brown fox jumps over the lazy dog"
strsplit(pangram, " ")
```

```
## [[1]]
## [1] "the"   "quick" "brown" "fox"   "jumps" "over"  "the"   "lazy"  "dog"
```

                
The output from strsplit is in a list. 
Notice that the unusual first line of `strsplit()`’s output consists of `[[1]]`. 
Similar to the way that R displays vectors, `[[1]]` means that R is showing the first element of a list. 
Lists are extremely important concepts in R; they allow you to combine all kinds of variables.
For example, a list can be made up of many elements, and elements could be vectors, data frames, matrices, or further lists. 

In this example, this list has only a single element. Yes, that’s right: The list has one element, but that element is a vector.

To extract an element from a list, you have to use double square brackets: `[[`. 
Split your pangram into words, and assign the first element to a new variable called words, using double-square-brackets `[[]]` subsetting, as follows:


```r
words <- strsplit(pangram, " ")[[1]]
```


We can then pull out the different words using indexing, where `words[1]` is the first element in the vector of `words` and `words[9]` would be the last:


```r
words[2]
```

```
## [1] "quick"
```

```r
words[9]
```

```
## [1] "dog"
```


1. Write a function `out()` that returns a string made up of just the first and last characters of its input.
    a. Outline the steps you need to take to write this function. Discuss with the person sitting next you.
    b. Write part of the code, make sure it works.
    c. Write the next step.
    d. Test your function.
    e. Can your function handle words of different lengths?


```r
out <- function(word) {
    letter <- strsplit(word, "")[[1]]
    abbrev <- paste(letter[1], letter[length(letter)], sep = "")
    abbrev
}

out("helium")
```

```
## [1] "hm"
```


__Making our function work with different inputs__. If we want just the last word, but we can't remember how long our sentence is, we can use `length()`


```r
length(words)
```

```
## [1] 9
```

```r
words[length(words)]
```

```
## [1] "dog"
```


__BREAK__

## Explaining the R Environments

Let's take a closer look at what happens when we call `fahr_to_celsius(32)`. To make things clearer, we'll start by putting the initial value 32 in a variable and store the final result in one as well:


```r
original <- 32
final <- fahr_to_celsius(original)
```


_Discuss and draw a diagram showing what memory looks like after the first line has been executed. Point to the environment_

When we call `fahr_to_celsius()`, R creates a new environment called the *evaluation environment* that is *local* to the function. This environment consists of two things

 1. a *frame*, which contains the symbol-value pairs (the assignment, called *binding*, of a value to a variable),
 2. an *enclosure*, which points to the enclosing environment, i.e. the environment where the function was called.

Initially, the *frame* only holds the object `temp`. Until is is used, `temp` in *unevaluated*; in effect it is a placeholder for whatever value we passed to the `temp` argument, which gets resolved as needed. This is known as *lazy evaluation*.

*Scoping* refers to the rules by which languages look up the value of a variable (symbol). R has two type of scoping rules

 1. *lexically* scoping, and
 2. *dynamic* scoping.

We won't discuss *dynamic* scoping heres. Basically, *lexical* scoping means that R looks up values for variables (symbols) based on how functions were nested when they were *called*. If a name isn't defined inside a function, R will look for the name in in the *parent frame*, the frame in the calling environment that is one where the function was called


```r
a <- 10
f <- function() {
    b <- 5
    a * b
}
f()
```

```
## [1] 50
```

There is no name `a` defined within the body of `f()`. Following lexical scoping rules, R will look in the *parent frame* of `f()` for a name `a`, which is found and has value `10`. As `f()` is defined in the workspace, the enclosing environment of `f()` is the *global environment*. Regardless of how functions are nested when called, R will always look for a name from the inside out until the global environment is reached.

When we nest a call to `fahr_to_kelvin()` inside `fahr_to_celsius()`, R creates another local environment to hold the variables involved in the evaluation of `fahr_to_kelvin()`.

When `fahr_to_celsius()` is called, it in turn calls `fahr_to_kelvin()`. As `temp` is passed to `fahr_to_kelvin()`, R evaluates `temp` to derive its value and this variable is passed to `fahr_to_kelvin()`. Now there are two `temp`s in play

 1. the `temp` local to `fahr_to_kelvin()`, and
 2. the `temp` local to `fahr_to_celcius()`.

As far as `fahr_to_kelvin()` is concerned there is only one `temp`; the `temp` local to it *masks* the other `temp`. That R knows which versions of variables with the same name belong to which function is due to environments.

When the call to `fahr_to_kelvin()` returns a value R throws away `fahr_to_kelvin()`'s frame and creates a new variable `temp_k` in the frame for `fahr_to_celsius()` to hold the temperature in Kelvin.

`kelvin_to_celsius()` is then called and a new environment to hold that function's variables is created.

Once again, R throws away that stack frame when `kelvin_to_celsius()` is done and creates the variable `result` in the frame of `fahr_to_celsius()`

Finally, when `fahr_to_celsius()` returns, R throws away its environment and puts `result` in a new variable called `final` that lives in the *global environment* where we called `fahr_to_celcius()` in the first place.

The summary of this is that the parent frame (the global environment is the parent frame of `fahr_to_celsius()`, the local environment of `fahr_to_celsius()` is the parent frame of `fahr_to_kelvin()`) is the environment where a function was defined (lexical scoping), the parent frame is the frame/environment from which the function was called.

The *global environment* is the final environment/frame and is always present; it holds the variables we defined outside the functions in our code. What it doesn't hold are the variables created during function calls. If we try to get the value of `temp` *after* our functions have finished running, R tells us that there's no such thing:


```r
paste("final value of temp after all function calls:", temp)
```

```
## Error: object 'temp' not found
```


Why go to all this trouble? Well, here's a function called `span()` that calculates the difference between the mininum and maximum values in an array:


```r
span <- function(a) {
    diff <- max(a) - min(a)
    diff
}
```


Notice `span()` assigns a value to variable called `diff`. We might very well use a variable with the same name (`diff`) to hold data:


```r
diff <- c(46, 55, 26, 64, 31, 68, 100, 79, 39, 95)
span(diff)
```

```
## [1] 74
```


We don't expect the variable `diff` to have the value 74 after this function call, so the name `diff` cannot refer to the same variable defined inside `span()` as it does in your workspace (the global environment). And yes, we could probably choose a different name than `diff` for our variable in this case, but we don't want to have to read every line of R code we write that we use to see what variable names its functions use before calling any of those functions, just in case they change the values of our variables.

The big idea here is __encapsulation__, and it's the key to writing correct, comprehensible programs. A function's job is to turn several operations into one so that we can think about a single function call instead of a dozen or a hundred statements each time we want to do something. That only works if functions don't interfere with each other; if they do, we have to pay attention to the details once again, which quickly overloads our short-term memory.

### Challenges

1. We previously wrote functions called `fence()` and `out()`. Walk your neighbor step-by-step through what happens when we call


```r
abbrev <- out(fence("carbon", "+"))
```


Use words or a diagram to explain how the variables change, and how the stack/environment changes. Look at the environment tab in RStudio.

__BREAK__

# Testing and Documenting

Once we start putting things in functions so that we can re-use them, we need to start testing that those functions are working correctly. To see how to do this, let's write a function to center a dataset around a particular value:


```r
center <- function(data, desired) {
    new <- (data - mean(data)) + desired
    new
}
```


We could test this on our actual data, but since we don't know what the values ought to be, it will be hard to tell if the result was correct. Instead, let's create a matrix of 0s and then center that around 3. This will make it simple to see if our function is working as expected:


```r
z <- matrix(nrow = 2, ncol = 2)
center(z, 3)
```

```
##      [,1] [,2]
## [1,]   NA   NA
## [2,]   NA   NA
```


That looks right, so let's try center on our real data:




```r
center(mat, 0)
```


It's hard to tell from the default output whether the result is correct, but there are a few simple tests that will reassure us:


```r
paste("original min, mean, and max are:", min(mat), ",", mean(mat), ",", max(mat))
```

```
## [1] "original min, mean, and max are: 1 , 800.5 , 1600"
```

```r
centered <- center(mat, 0)

paste("original min, mean, and max are:", min(centered), ",", mean(centered), 
    ",", max(centered))
```

```
## [1] "original min, mean, and max are: -799.5 , 0 , 799.5"
```


That seems almost right: the original mean was about 6.1, so the lower bound from zero is how about -6.1. The mean of the centered data isn't quite zero --- we'll explore why not in the challenges --- but it's pretty close. We can even go further and check that the standard deviation hasn't changed:


```r
paste("std dev before and after:", sd(mat), sd(centered))
```

```
## [1] "std dev before and after: 462.024530373298 462.024530373298"
```


Those values look the same, but we probably wouldn't notice if they were different in the sixth decimal place. Let's do this instead:


```r
paste("difference in standard deviations before and after:", sd(mat) - sd(centered))
```

```
## [1] "difference in standard deviations before and after: 0"
```


Sometimes, a very small difference can be detected. This could be due to rounding at very low decimal places. R has a useful function for comparing two objects allowing for rounding errors, `all.equal()`:


```r
all.equal(sd(mat), sd(centered))
```

```
## [1] TRUE
```


It's still possible that our function is wrong, but it seems unlikely enough that we should probably get back to doing our analysis. We have one more task first, though: we should write some documentation for our function to remind ourselves later what it's for and how to use it.

A common way to put documentation in software is to add comments like this:


```r
# return a new matrix containing the original data centered around the
# desired value.
center <- function(data, desired) {
    new <- (data - mean(data)) + desired
    new
}
```


Formal documentation for R functions is written in separate `.Rd` using a markup language similar to LaTeX. The **roxygen2** package allows R coders to write documentation alongside the function code and then process it into the appropriate `.Rd` files.

### Challenges

This next challenge has several steps. Think about how you break down a difficult problem into manageable pieces.

1. Write a function called `analyze()` that takes a filename as a parameter and displays the 3 graphs you made earlier (average, min and max inflammation over time). i.e., `analyze("data/inflammation-01.csv")` should produce the graphs already shown, while `analyze("inflammation-02.csv")` should produce corresponding graphs for the second data set. Be sure to give your function a docstring.

## Defining Defaults

We have passed parameters to functions in two ways: directly, as in `dim(mat)`, and by name, as in `matrix(data = 0, nrow = 2, ncol = 2)`. We can pass arguments to functions without naming them


```r
matrix(0, 2, 2)
```

```
##      [,1] [,2]
## [1,]    0    0
## [2,]    0    0
```


To understand what's going on, and make our own functions easier to use, let's re-define our `center()` function like this:


```r
center <- function(data, desired = 0) {
    # return a new matrix containing the original data centered around the
    # desired value.
    new <- (data - mean(data)) + desired
    new
}
```


The key change is that the second argument is now written `desired = 0` instead of just `desired`. If we call the function with two arguments, it works as it did before:


```r
test_data <- matrix(0, 2, 2)
center(test_data, 3)
```

```
##      [,1] [,2]
## [1,]    3    3
## [2,]    3    3
```


But we can also now call `center()` with just one parameter, in which case `desired` is automatically assigned the default value of `0`:


```r
more_data <- matrix(0, 2, 2) + 5
more_data
```

```
##      [,1] [,2]
## [1,]    5    5
## [2,]    5    5
```

```r
center(more_data)
```

```
##      [,1] [,2]
## [1,]    0    0
## [2,]    0    0
```


This is handy: if we usually want a function to work one way, but occasionally need it to do something else, we can allow people to pass a parameter when they need to but provide a default to make the normal case easier.

R has three ways that arguments supplied by you are matched to the *formal arguments * of the function definition

 1. by complete name, 
 2. by partial name (matching on initial *n* characters of the argument name), and
 3. by position.

Arguments are matched in the manner outlined above in *that order*, by complete name, then by partial matching of names, and finally by position.

The example below shows how R matches values to arguments


```r
display <- function(a = 1, b = 2, c = 3) {
    paste("a:", a, "b:", b, "c:", c)
}

paste("no parameters:", display())
```

```
## [1] "no parameters: a: 1 b: 2 c: 3"
```

```r
paste("one parameter:", display(55))
```

```
## [1] "one parameter: a: 55 b: 2 c: 3"
```

```r
paste("two parameters:", display(55, 66))
```

```
## [1] "two parameters: a: 55 b: 66 c: 3"
```

```r
paste("three parameters:", display(55, 66, 77))
```

```
## [1] "three parameters: a: 55 b: 66 c: 77"
```


As this example shows, arguments are matched from left to right, and any that haven't been given a value explicitly get their default value. We can override this behavior by naming the value as we pass it in:


```r
paste("only setting the value of c", display(c = 77))
```

```
## [1] "only setting the value of c a: 1 b: 2 c: 77"
```


With that in hand, let's look at the help for `read.csv()`:


```r
`?`(read.csv)
```


_walk through the help file, point out important items. How to read the defaults, how to know what you need to specify. Definitions and examples..._

There's a lot of information here, but the most important part is the first couple of lines:


```r
read.csv(file, header = TRUE, sep = ",", quote = "\"", dec = ".", fill = TRUE, 
    comment.char = "", ...)
```


This tells us that `read.csv()` has one argument, `file`, that doesn't have a default value, and eight others that do. If we call the function like this:


```r
read.csv("data/inflammation-01.csv", ",")
```

```
## Error: invalid argument type
```


the filename is assigned to`file` (which is what we want), but the delimiter string `","` is assigned to the argument `header` rather than `sep`, because `header` is the second parameter in the list. That's why we don't have to provide `file =` for the filename, but do have to provide `sep =` for the second parameter.


```r
read.csv("data/inflammation-01.csv", sep = ",")
```


### Challenges

Rewrite the `center()` function so that it scales data to lie between 0.0 and 1.0 by default, but will allow the caller to specify lower and upper bounds if they want. Compare your implementation to your neighbor's: do the two functions always behave the same way?

Key Points
===============
* Define a function using `function` name(...params...).
* The body of a function should be indented.
* Call a function using name(...values...).
* Numbers are stored as integers or floating-point numbers.
* Each time a function is called, a new stack frame is created on the call stack to hold its parameters and local variables.
* R looks for variables in the current environment before looking for them at the top level.
* Use help(thing) to view help for something.
* Put docstrings in functions to provide help for that function.
* Annotate your code!
* Specify default values for parameters when defining a function using name=value in the parameter list.
* Parameters can be passed by matching based on name, by position, or by omitting them (in which case the default value is used).

## Next Steps

We now have a function called analyze to visualize a single data set. We could use it to explore all 12 of our current data sets like this:


```r
analyze("data/inflammation-01.csv")
analyze("data/inflammation-02.csv")
# ...
analyze("data/inflammation-12.csv")
```


but the chances of us typing all 12 filenames correctly aren't great, and we'll be even worse off if we get another hundred files. What we need is a way to tell R to do something once for each file, and that will be the subject of the next lesson.
