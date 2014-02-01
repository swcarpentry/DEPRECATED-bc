---
title: R Objects
---

All R code deals with objects. An object is a *thing* that is
represented by the computer. Data can be represented by different
objects of a certain type and mode.

The *type* or *mode* of an object defines how it is stored. It could
be

- a character,
- a numeric value,
- an integer,
- a complex number, or
- a logical value (Boolean value: TRUE/FALSE).

It can be determined with the function `mode()`.

Objects also have a *class* or *data structure*: what information is
hold by the object and how the object may be used. It could be

- a vector,
- a list,
- a matrix,
- an array, or 
- a data frame.

It can be determined with the function `class()`.

## Vectors
It is the most basic data type. Vectors only have one dimension and
all their elements must be the same mode. There are various ways to
create vectors. The simplest one is with the `c` function.

    v1 <- c(1, 2, 5)

The `c` function *coerces* all of its argument into a single type:

    v2 <- c(1, 3, "a")
    mode(v2)

Objects can be explicitly coerced with the `as.` function.

    as.character(v1)

You can also use the `:` operator or the `seq` function:

    1:10
    seq(from = 5, to = 25, by = 5)

### Factors
Factor is a special type of vector to store categorical values.

    breed <- factor(c("Holstein", "Brown Swiss", "Holstein",
	"Ayrshire", "Canadian"))

It stores values as a set of labeled integers. Some functions treat
factors differently from numeric vectors.

    table(breed)

## Lists
A list is an ordered collection of objects where the objects can be of
different modes.

    l <- list("a", "b", "c")

Each element of a list can be given a name and referred to by that
name. Elements of a list can be accessed by their number or their name.

    cow <- list(breed = "Holstein", age = 3, last.prod = c(25, 35, 32))
    cow$breed
	cow[[1]]

Lists can be used to hold together multiple values returned from a
function. For example the elements used to create an histogram can be
saved and returned:

    h <- hist(islands)
	str(h)

The function `str` is used here. It stands for *structure* and shows
the internal structure of an R object.

## Matrices
A matrix is a 2-dimensional vector, of a single mode. You can create a
matrix with the `matrix` function:

	(m <- rbind(c(1, 4), c(2, 2)))
	m <- matrix(data = 1:12,
		nrow = 4, ncol = 3,
		dimnames = list(c("cow1", "cow2", "cow3", "cow4"),
			c("milk", "fat", "prot")))

## Arrays
An array is a multidimensional vector, of a single type. It can be
created with the `array` function:

    a <- array(data = 1:24, dim = c(3, 4, 2))

## Data frames
Data frames are used to store tabular data: multiple rows, columns and
format.

    df <- data.frame(cow = c("Moo-Moo", "Daisy", "Elsie"),
		prod = c(35, 40, 28),
		pregnant = c(TRUE, FALSE, TRUE))

