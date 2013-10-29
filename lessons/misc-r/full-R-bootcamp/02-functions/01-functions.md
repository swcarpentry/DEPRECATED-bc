
# Writing functions in R

If you have to repeat the same few lines of code more than once, then you really need to write a function. Functions are a fundamental building block of R. You use them all the time in R and it's not that much harder to string functions together (or write entirely new ones from scratch) to do more.

* R functions are objects just like anything else. 
* By default, R function arguments are lazy - they're only evaluated if they're actually used:
* Every call on a R object is almost always a function call.

## Basic components of a function

* The `body()`, the code inside the function.
* The `formals()`, the "formal" argument list, which controls how you can call the function.
* The `environment()`` which determines how variables referred to inside the function are found.
* `args()` to list arguments.

```
f <- function(x) x
f

formals(f)

environment(f)
```

**Question: How do we delete this function from our environment?**

## More on environments
Variables defined inside functions exist in a different environment than the global environment. However, if a function is not defined inside one, it will look one level above.

example.

```
x <- 2
g <- function() { 
  y <- 1
  c(x, y)
}  
g()
rm(x, g)
```

Same rule applies for nested functions.




A first useful function.

```coffee
first <- function(x, y) {
    z <- x + y
    return(z)
}
```

```
add <- function(a, b) {
  return (a + b)
}
vector <- c(3, 4, 5, 6)

sapply(vector, add, 1)
```

## What does this function return?

```coffee
x <- 5
f <- function() {
y <- 10
c(x = x, y = y) }
f()
```

## What does this function return?

```coffee
x <- 5
g <- function() {
  x <- 20
  y <- 10
  c(x = x, y = y)
} 
g()
```

## What does this function return??

```coffee
x <- 5
h <- function() {
  y <- 10
  i <- function() {
z <- 20
    c(x = x, y = y, z = z)
  }
i() 
}
h()
```


## Functions with pre defined values

```coffee
temp <- function(a = 1, b =2) {
    return(a + b)
}
```

## Functions usually return the last value it computed

```
f <- function(x) {
  if (x < 10) {
    0
  } else {
    10
  }
}
f(5)
f(15)
```

