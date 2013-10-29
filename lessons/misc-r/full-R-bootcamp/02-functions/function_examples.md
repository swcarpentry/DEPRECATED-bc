# Functions 

Represent one of the most powerful features of the R language. By creating functions, you can go past just using what others have written and create your own.

After a while operations at the interactive prompt are not enough. You have special operations that need to be grouped together. That's when you'll need to write functions.

* Basic of how to write functions
* Lexical scoping and scoping rules
* Example

Functions are created using the `function()` directive.
They are stored as any other R object.  They are of class "function"

**Basic syntax**

```
func <- function(<take some input>)  {
    # run a whole bunch of your code
}
```

Functions are first class objects. Treat them the same way as any other objects. This means you can pass them to other functions (as you would do with say vectors). Very useful feature for statistical analyses. We can also nest functions inside each other.

Return value is the very last R expression to be evaluated. No special expression for returning although you can use return()

Functions have named arguments and these can have default values. Arguments used inside a function definition are called formal arguments.

use `formals()` to see what they are.

Not every function makes use of all formal arguments.

Function arguments can be missing or might just use default values. Formal arguments are included in the function definition.

If a function has 10 arguments, you don't have to specify all of them. 

R function arguments can be matched positionally or by name.

```
args(sd)
```
Takes x, a vector of data. 
Default is `na.rm` is false.
So you don't have to do anything.

```
dat <- rnorm(1000)
sd(dat) # Matched positionally
sd(x = dat) # matched by name
sd(x = dat, na.rm = FALSE)
# When naming you can switch them but not recommended
sd(na.rm = FALSE, dat)
```

Even though it's allowed, don't switch positions on names.

See 

```
args(lm)
```

It's good to name arguments when you have a long list to work through. Arguments can also be partially matched.

```
add <- function(x = 1, y =2) {
x + y
}
```

R does what's known as lazy evaluation. Arguments to functions are lazily evaluated. Common model in many languages. Only evaluated as they are needed.

```
add <- function(a, b) {
a^2
}
```

The function call never uses b. So calling `f(10)` will not produce an error because `a` got positionally matched.

Another example of lazy evaluation


```
add <- function(a, b) {
print(a^2)
print(b^2)
}
```

`add(10)` will work until it reaches a point where the function will break.

**The three dot operator.**

`...` are used when extending another function and you don't want to copy the entire list of arguments from the original function.

Also useful when the number of arguments cannot be known in advance.

e.g. `paste()`

takes variable number of arguments. Function does not know how many things it's going to paste together.

