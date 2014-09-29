---
layout: lesson
root: ../..
title: R Reference
---
# Basic Operation

- `# this is a comment in R`
- Use `x <- 3` to assign a value, `3`,  to a variable, `x`
- R counts from 1, unlike many other programming languages (e.g., Python)
- `length(thing)` produces the length of a collection
- `c(value1, value2, value3)` creates a vector
- `container[i]` selects the i'th element from the container

List objects in current environment
`ls()`

Remove objects in current environment
`rm(x)`

Remove all objects from current environment
`rm(list = ls())`

# Control Flow

- Create a contitional using `if`, `elif`, and `else`

		if(x > 0){
			print("value is positive")
		} else if (x < 0){
			print("value is negative")
		} else{
			print("value is neither positive nor negative")
		}

- create a `for` loop to process elements in a collection one at a time

		for (i in 1:5) {
			print(i)
		}

This will print:

		1
		2
		3
		4
		5


- Use `==` to test for equality
  - `3 == 3`, will return `TRUE`,
  - `'apple' == 'orange'` will return `FALSE`
- `X & Y` is `True` is both X and Y are true
- `X | Y` is `True` if either X or Y, or both are true

# Functions

- Defining a function:

		is_positive <- function(integer_value){
			if(interver_value > 0){
			   TRUE
			else{
			   FALSE
			{
		}

In R, the last executed line of a function is automatically returned

- Specifying a default value for a function argrment

		increment_me <- function(value_to_increment, value_to_increment_by = 1){
			value_to_increment + value_to_increment_by
		}

`increment_me(4)`, will return 5

`intrement_me(4, 6)`, will return 10

- Call a function by using `function_name(function_arguments)`

	- apply family of functions:

			apply()
			sapply()
			lapply()
			mapply()

`apply(dat, MARGIN = 2, mean)`
will return the average (`mean`) of each column in `dat`

# Packages
- Install package by using `install.packages("package-name")`
- Update packages by using `update.packages("package-name")`
- Load packages by using `library("package-name")`
