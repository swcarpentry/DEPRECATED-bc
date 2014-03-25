Novice R materials - following Python examples
--------------------------------------------------

* Draw concept map of where we are headed towards best scientific practices, and reproducibility.
* Its really important that you know what you did. More journals/grants/etc. are also making it important for them to know what you did.
* A lot of scientific code is NOT reproducible.
* If you keep a lab notebook, why are we not as careful with our code. 
* We edit each others manuscripts, but we don't edit each other's code. 
* If you write your code with "future you" in mind, you will save yourself and others a lot of time.

Very basics of R
-----------------
R is a versatile, open source programming/scripting language that's useful both for statistics but also data science. Inspired by the programming language S.

* Open source software under GPL.
* Superior (if not just comparable) to commercial alternatives. R has over 5,000 user contributed packages at this time. It's widely used both in academia and industry.
* Available on all platforms.
* Not just for statistics, but also general purpose programming.
* Is object oriented and functional.
* Large and growing community of peers.

Some of the same commands we learned from the command line can be used in R.
List objects in your current environment

	ls()

Remove objects from your current environment.

	x <- 5
	rm(x)
	x
	
Remove all objects from your current environment. 

	rm(list = ls())
Notice that we have _nested_ one function inside of another.

Use # signs to comment. Comment liberally in your R scripts. Anything to the right of a # is ignored by R.

__Assignment operator__

`<-` is the assignment operator. Assigns values on the right to objects on the left. Mostly similar to `=` but not always. Learn to use <- as it is good programming practice. Using = in place of <- can lead to issues down the line.

__Package management__

`install.packages("package-name")` will download a package from one of the CRAN mirrors assuming that a binary is available for your operating system. If you have not set a preferred CRAN mirror in your options(), then a menu will pop up asking you to choose a location.

Use `old.packages()` to list all your locally installed packages that are now out of date. update.packages() - will update all packages in the known libraries interactively. This can take a while if you haven't done it recently. To update everything without any user intervention, use the ask = FALSE argument.

	update.packages(ask = FALSE)



Introduction to R and RStudio
--------------------------------
Let's start by learning about our tool. 

_Point out the different windows in R._ 
* Console, Scripts, Environments, Plots
* Avoid using shortcuts. 
* Code and workflow is more reproducible if we can document everything that we do.
* Our end goal is not just to "do stuff" but to do it in a way that anyone can easily and exactly replicate our workflow and results.

You can get output from R simply by typing in math
	
	3 + 5
	12/7

or by typing words, with the command "print"
	
	print ("hello world")

We can annotate our code (take notes) by typing "#". Everything to the right of # is ignored by R

We can save our results to an object, if we give it a name
	
	a = 60 * 60
	hours <- 365 * 24

Data types and structures
----------------------------
__Understanding basic data types in R__

To make the best of the R language, you'll need a strong understanding of the basic data types and data structures and how to operate on those.

Very Important to understand because these are the objects you will manipulate on a day-to-day basis in R. Dealing with object conversions is one of the most common sources of frustration for beginners.

Everything in R is an object.

R has 6 (although we will not discuss the raw class for this workshop) atomic classes.

* character
* numeric (real or decimal)
* integer
* logical
* complex

__Example	Type__

* “a”, “swc”	character
* 2, 15.5	numeric
* 2 (Must add a L at end to denote integer)	integer
* TRUE, FALSE	logical
* 1+4i	complex


		typeof() # what is it?
		length() # how long is it? What about two dimensional objects?
		attributes() # does it have any metadata?

		# Example
		x <- "dataset"
		typeof(x)
		attributes(x)

		y <- 1:10
		typeof(y)
		length(y)
		attributes(y)

		z <- c(1L, 2L, 3L)
		typeof(z)

R has many __data structures__. These include

* atomic vector
* list
* matrix
* data frame
* factors
* tables


Vectors
--------

A vector is the most common and basic data structure in `R` and is pretty much the workhorse of R. Technically, vectors can be one of two types:

* atomic vectors
* lists

although the term "vector" most commonly refers to the atomic type not lists.


**Atomic Vectors**

A vector can be a vector of elements that are most commonly `character`, `logical`, `integer` or `numeric`.

You can create an empty vector with `vector()` (By default the mode is `logical`. You can be more explicit as shown in the examples below.) It is more common to use direct constructors such as `character()`, `numeric()`, etc.


	x <- vector()
	# with a length and type
	vector("character", length = 10)
	character(5) ## character vector of length 5
	numeric(5)
	logical(5)

Various examples:

	x <- c(1, 2, 3)
	x
	length(x)

`x` is a numeric vector. These are the most common kind. They are numeric objects and are treated as double precision real numbers. To explicitly create integers, add an `L` at the end.

	x1 <- c(1L, 2L, 3L)

You can also have logical vectors. 

	y <- c(TRUE, TRUE, FALSE, FALSE)

Finally you can have character vectors:

	z <- c("Sarah", "Tracy", "Jon")

**Examine your vector**  

	typeof(z)
	length(z)
	class(z)
	str(z)

Question: Do you see a property that's common to all these vectors above?

**Add elements**

	z <- c(z, "Annette")
	z


More examples of vectors

	x <- c(0.5, 0.7)
	x <- c(TRUE, FALSE)
	x <- c("a", "b", "c", "d", "e")
	x <- 9:100
	x <- c(1+0i, 2+4i)

You can also create vectors as a sequence of numbers

	series <- 1:10
	seq(10)
	seq(1, 10, by = 0.1)
	
`Inf` is infinity. You can have either positive or negative infinity.

	1/0

`NaN` means Not a number. It's an undefined value.
	
	0/0

Each object can have __attributes__. Attribues can be part of an object of R. These include:

* names
* dimnames
* dim
* class
* attributes (contain metadata)

You can also glean other attribute-like information such as length (works on vectors and lists) or number of characters (for character strings).

	length(1:10)

	nchar("Software Carpentry")

What happens when you mix types?

R will create a resulting vector that is the least common denominator. The coercion will move towards the one that's easiest to __coerce__ to.

Guess what the following do without running them first

	xx <- c(1.7, "a") 
	xx <- c(TRUE, 2) 
	xx <- c("a", TRUE) 
	
This is called implicit coercion. You can also coerce vectors explicitly using the `as.<class_name>`. Example

	as.numeric()
	as.character()

Matrix
---------

Matrices are a special vector in R. They are not a separate type of object but simply an atomic vector with dimensions added on to it. Matrices have rows and columns.

	m <- matrix(nrow = 2, ncol = 2)
	m

	dim(m)

Matrices are filled column-wise.

	m <- matrix(1:6, nrow = 2, ncol = 3)

Other ways to construct a matrix

	m <- 1:10
	dim(m) <- c(2, 5)

This takes a vector and transform into a matrix with 2 rows and 5 columns.

Another way is to bind columns or rows using cbind() and rbind().

	x <- 1:3
	y <- 10:12
	cbind(x, y)

	rbind(x, y)

You can also use the byrow argument to specify how the matrix is filled. From R's own documentation:

	mdat <- matrix(c(1,2,3, 11,12,13), nrow = 2, ncol = 3, byrow = TRUE,
               dimnames = list(c("row1", "row2"),
                               c("C.1", "C.2", "C.3")))
	mdat
	
List
-----

In R lists act as containers. Unlike atomic vectors, the contents of a list are not restricted to a single mode and can encompass any mixture of data types. Lists are sometimes called recursive vectors, because a list can contain other lists. This makes them fundamentally different from atomic vectors.

A list is a special type of vector. Each element can be a different type.

Create lists using `list()` or coerce other objects using `as.list()`

	x <- list(1, "a", TRUE, 1+4i)
	x

	x <- 1:10
	x <- as.list(x)
	length(x)

1. What is the class of x[1]?
2. How about x[[1]]?

	xlist <- list(a = "Karthik Ram", b = 1:10, data = head(iris))
	xlist

1. What is the length of this object? What about its structure?

Lists can be extremely useful inside functions. You can “staple” together lots of different kinds of results into a single object that a function can return.

A list does not print to the console like a vector. Instead, each element of the list starts on a new line.

Elements are indexed by double brackets. Single brackets will still return a(nother) list.


Factors
----------

Factors are special vectors that represent categorical data. Factors can be ordered or unordered and are important for modelling functions such as lm() and glm() and also in plot methods.

Factors can only contain pre-defined values.

Factors are pretty much integers that have labels on them. While factors look (and often behave) like character vectors, they are actually integers under the hood, and you need to be careful when treating them like strings. Some string methods will coerce factors to strings, while others will throw an error.

Sometimes factors can be left unordered. Example: male, female.

Other times you might want factors to be ordered (or ranked). Example: low, medium, high.

Underlying it's represented by numbers 1, 2, 3.

They are better than using simple integer labels because factors are what are called self describing. male and female is more descriptive than 1s and 2s. Helpful when there is no additional metadata.

Which is male? 1 or 2? You wouldn't be able to tell with just integer data. Factors have this information built in.

Factors can be created with `factor()`. Input is generally a character vector.

	x <- factor(c("yes", "no", "no", "yes", "yes"))
	x

`table(x)` will return a frequency table.

If you need to convert a factor to a character vector, simply use

	as.character(x)

In modeling functions, it is important to know what the baseline level is. This is the first factor but by default the ordering is determined by alphabetical order of words entered. You can change this by speciying the levels (another option is to use the function relevel).

	x <- factor(c("yes", "no", "yes"), levels = c("yes", "no"))
	x

Data frame
------------

A data frame is a very important data type in R. It's pretty much the de facto data structure for most tabular data and what we use for statistics.

Data frames can have additional attributes such as rownames(), which can be useful for annotating data, like subject_id or sample_id. But most of the time they are not used.

Some additional information on data frames:

* Usually created by read.csv() and read.table().
* Can convert to matrix with data.matrix()
* Coercion will be forced and not always what you expect.
* Can also create with data.frame() function.
* Find the number of rows and columns with nrow(df) and ncol(df), respectively.
* Rownames are usually 1..n.

__Combining data frames__
-------------

	df <- data.frame(id = letters[1:10], x = 1:10, y = 11:20)
	df

__Useful functions__
-------------

* head() - see first 6 rows
* tail() - see last 6 rows
* dim() - see dimensions
* nrow() - number of rows
* ncol() - number of columns
* str() - structure of each column
* names() - will list the names attribute for a data frame (or any object really), which gives the column names.
* A data frame is a special type of list where every element of the list has same length.

See that it is actually a special list:

	is.list(iris)
	class(iris)


| Dimensions | Homogenous | Heterogeneous |
| ------- | ---- | ---- |
| 1-D | atomic vector | list |
| 2_D | matrix | dataframe |



__Indexing__
-----------

Vectors have positions, these positions are ordered and can be called using name_vector[index]
	
	names[2]

__Functions__
---------

A function is a saved object that takes inputs to perform a task. 
Functions take in information and return desired outputs.

output = name_of_function(inputs)
	
	y = sum(x)

__Help__
--------

All functions come with a help screen. 
It is critical that you learn to read the help screens since they provide important information on what the function does, 
how it works, and usually sample examples at the very bottom.

__Install new functions__
--------

To install any new package install.packages('ggplot2')

You can't ever learn all of R, but you can learn how to build a program and how to find help
to do the things that you want to do. Let's get hands-on.