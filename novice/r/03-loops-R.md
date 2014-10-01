---
layout: lesson
root: ../..
---



## Analyzing Multiple Data Sets

We have created a function called `analyze` that creates graphs of the minimum, average, and maximum daily inflammation rates for a single data set:


<pre class='in'><code>analyze <- function(filename) {
  # Plots the average, min, and max inflammation over time.
  # Input is character string of a csv file.
  dat <- read.csv(file = filename, header = FALSE)
  avg_day_inflammation <- apply(dat, 2, mean)
  plot(avg_day_inflammation)
  max_day_inflammation <- apply(dat, 2, max)
  plot(max_day_inflammation)
  min_day_inflammation <- apply(dat, 2, min)
  plot(min_day_inflammation)
}

analyze("inflammation-01.csv")</code></pre>

<img src="figure/03-loops-R-inflammation-011.png" title="plot of chunk inflammation-01" alt="plot of chunk inflammation-01" style="display: block; margin: auto;" /><img src="figure/03-loops-R-inflammation-012.png" title="plot of chunk inflammation-01" alt="plot of chunk inflammation-01" style="display: block; margin: auto;" /><img src="figure/03-loops-R-inflammation-013.png" title="plot of chunk inflammation-01" alt="plot of chunk inflammation-01" style="display: block; margin: auto;" />

We can use it to analyze other data sets one by one:


<pre class='in'><code>analyze("inflammation-02.csv")</code></pre>

<img src="figure/03-loops-R-inflammation-021.png" title="plot of chunk inflammation-02" alt="plot of chunk inflammation-02" style="display: block; margin: auto;" /><img src="figure/03-loops-R-inflammation-022.png" title="plot of chunk inflammation-02" alt="plot of chunk inflammation-02" style="display: block; margin: auto;" /><img src="figure/03-loops-R-inflammation-023.png" title="plot of chunk inflammation-02" alt="plot of chunk inflammation-02" style="display: block; margin: auto;" />

but we have a dozen data sets right now and more on the way.
We want to create plots for all our data sets with a single statement.
To do that, we'll have to teach the computer how to repeat things.

<div class="objectives" markdown="1">
#### Objectives

* Explain what a `for` loop does.
* Correctly write `for` loops to repeat simple calculations.
* Trace changes to a loop variable as the loop runs.
* Trace changes to other variables as they are updated by a `for` loop.
* Use a function to get a list of filenames that match a simple pattern.
* Use a `for` loop to process multiple files.
</div>

### For Loops

Suppose we want to print each word in a sentence.
One way is to use six `print` statements:


<pre class='in'><code>best_practice <- c("Let", "the", "computer", "do", "the", "work")
print_words <- function(sentence) {
  print(sentence[1])
  print(sentence[2])
  print(sentence[3])
  print(sentence[4])
  print(sentence[5])
  print(sentence[6])
}

print_words(best_practice)</code></pre>



<div class='out'><pre class='out'><code>[1] "Let"
[1] "the"
[1] "computer"
[1] "do"
[1] "the"
[1] "work"
</code></pre></div>

but that's a bad approach for two reasons:

 1. It doesn't scale: if we want to print the elements in a vector that's hundreds long, we'd be better off just typing them in.

 2. It's fragile: if we give it a longer vector, it only prints part of the data, and if we give it a shorter input, it returns `NA` values because we're asking for elements that don't exist!


<pre class='in'><code>best_practice[-6]</code></pre>



<div class='out'><pre class='out'><code>[1] "Let"      "the"      "computer" "do"       "the"     
</code></pre></div>



<pre class='in'><code>print_words(best_practice[-6])</code></pre>



<div class='out'><pre class='out'><code>[1] "Let"
[1] "the"
[1] "computer"
[1] "do"
[1] "the"
[1] NA
</code></pre></div>

> **Tip:** R has has a special variable, `NA`, for designating missing values that are **N**ot **A**vailable in a data set.
See `?NA` and [An Introduction to R][na] for more details.

[na]: http://cran.r-project.org/doc/manuals/r-release/R-intro.html#Missing-values

Here's a better approach:


<pre class='in'><code>print_words <- function(sentence) {
  for (word in sentence) {
    print(word)
  }
}

print_words(best_practice)</code></pre>



<div class='out'><pre class='out'><code>[1] "Let"
[1] "the"
[1] "computer"
[1] "do"
[1] "the"
[1] "work"
</code></pre></div>

This is shorter---certainly shorter than something that prints every character in a hundred-letter string---and more robust as well:


<pre class='in'><code>print_words(best_practice[-6])</code></pre>



<div class='out'><pre class='out'><code>[1] "Let"
[1] "the"
[1] "computer"
[1] "do"
[1] "the"
</code></pre></div>

The improved version of `print_words` uses a [for loop](../../gloss.html#for-loop) to repeat an operation---in this case, printing---once for each thing in a collection.
The general form of a loop is:


<pre class='in'><code>for (variable in collection) {
  do things with variable
}</code></pre>

We can name the [loop variable](../../gloss.html#loop-variable) anything we like (with a few [restrictions][], e.g. the name of the variable cannot start with a digit).
`in` is part of the `for` syntax.
Note that the body of the loop is enclosed in curly braces `{ }`.
For a single-line loop body, as here, the braces aren't needed, but it is good practice to include them as we did.

[restrictions]: http://cran.r-project.org/doc/manuals/R-intro.html#R-commands_003b-case-sensitivity-etc

Here's another loop that repeatedly updates a variable:


<pre class='in'><code>len <- 0
vowels <- c("a", "e", "i", "o", "u")
for (v in vowels) {
  len <- len + 1
}
# Number of vowels
len</code></pre>



<div class='out'><pre class='out'><code>[1] 5
</code></pre></div>

It's worth tracing the execution of this little program step by step.
Since there are five elements in the vector `vowels`, the statement inside the loop will be executed five times.
The first time around, `len` is zero (the value assigned to it on line 1) and `v` is `"a"`.
The statement adds 1 to the old value of `len`, producing 1, and updates `len` to refer to that new value.
The next time around, `v` is `"e"` and `len` is 1, so `len` is updated to be 2.
After three more updates, `len` is 5; since there is nothing left in the vector `vowels` for R to process, the loop finishes.

Note that a loop variable is just a variable that's being used to record progress in a loop.
It still exists after the loop is over, and we can re-use variables previously defined as loop variables as well:


<pre class='in'><code>letter <- "z"
for (letter in c("a", "b", "c")) {
  print(letter)
}</code></pre>



<div class='out'><pre class='out'><code>[1] "a"
[1] "b"
[1] "c"
</code></pre></div>



<pre class='in'><code># after the loop, letter is
letter</code></pre>



<div class='out'><pre class='out'><code>[1] "c"
</code></pre></div>

Note also that finding the length of a vector is such a common operation that R actually has a built-in function to do it called `length`:


<pre class='in'><code>length(vowels)</code></pre>



<div class='out'><pre class='out'><code>[1] 5
</code></pre></div>

`length` is much faster than any R function we could write ourselves, and much easier to read than a two-line loop; it will also give us the length of many other things that we haven't met yet, so we should always use it when we can (see this [lesson](00-first-timers.html) to learn more about the different ways to store data in R).

#### Challenges

1. R has a built-in function called `seq` that creates a list of numbers:


<pre class='in'><code>seq(3)</code></pre>



<div class='out'><pre class='out'><code>[1] 1 2 3
</code></pre></div>

Using `seq`, write a function that prints the first **N** natural numbers, one per line:




<pre class='in'><code>print_N(3)</code></pre>



<div class='out'><pre class='out'><code>[1] 1
[1] 2
[1] 3
</code></pre></div>

2. Exponentiation is built into R:


<pre class='in'><code>2^4</code></pre>



<div class='out'><pre class='out'><code>[1] 16
</code></pre></div>

Write a function called `expo` that uses a loop to calculate the same result.




<pre class='in'><code>expo(2, 4)</code></pre>



<div class='out'><pre class='out'><code>[1] 16
</code></pre></div>

1. Write a function called `total` that calculates the sum of the values in a vector.
(R has a built-in function called `sum` that does this for you.
Please don't use it for this exercise.)




<pre class='in'><code>ex_vec <- c(4, 8, 15, 16, 23, 42)
total(ex_vec)</code></pre>



<div class='out'><pre class='out'><code>[1] 108
</code></pre></div>

### Processing Multiple Files

We now have almost everything we need to process all our data files. 
The only thing that's missing is a function that finds files whose names match a pattern.
We do not need to write it ourselves because R already has a function to do this called `list.files`.

If we run the function without any arguments, `list.files()`, it returns every file in the current working directory.
We can understand this result by reading the help file (`?list.files`).
The first argument, `path`, is the path to the directory to be searched, and it has the default value of `"."` (recall from the [lesson](../shell/01-filedir.html) on the Unix Shell that `"."` is shorthand for the current working directory).
The second argument, `pattern`, is the pattern being searched, and it has the default value of `NULL`.
Since no pattern is specified to filter the files, all files are returned.

So to list all the csv files, we could run either of the following:


<pre class='in'><code>list.files(pattern = "csv")</code></pre>



<div class='out'><pre class='out'><code> [1] "inflammation-01.csv" "inflammation-02.csv" "inflammation-03.csv"
 [4] "inflammation-04.csv" "inflammation-05.csv" "inflammation-06.csv"
 [7] "inflammation-07.csv" "inflammation-08.csv" "inflammation-09.csv"
[10] "inflammation-10.csv" "inflammation-11.csv" "inflammation-12.csv"
</code></pre></div>



<pre class='in'><code>list.files(pattern = "inflammation")</code></pre>



<div class='out'><pre class='out'><code> [1] "inflammation-01.csv" "inflammation-02.csv" "inflammation-03.csv"
 [4] "inflammation-04.csv" "inflammation-05.csv" "inflammation-06.csv"
 [7] "inflammation-07.csv" "inflammation-08.csv" "inflammation-09.csv"
[10] "inflammation-10.csv" "inflammation-11.csv" "inflammation-12.csv"
</code></pre></div>

We have to name the argument `pattern` because we are using the default option for the first argument `path` (see the previous [lesson](02-func-R.html) to review default function arguments).

> **Tip:** Since this is just a small example, we have the data and code in the same directory.
For larger projects, it is recommended to organize separate parts of the analysis into multiple subdirectories, e.g. one subdirectory for the raw data, one for the code, and one for the results like figures.
For more advice on this topic, you can read [A quick guide to organizing computational biology projects][Noble2009] by William Stafford Noble.

[Noble2009]: http://www.ploscompbiol.org/article/info%3Adoi%2F10.1371%2Fjournal.pcbi.1000424

As these examples show, `list.files` result is a vector of strings, which means we can loop over it to do something with each filename in turn.
In our case, the "something" we want is our `analyze` function.
Let's test it by analyzing the first three files in the vector:


<pre class='in'><code>filenames <- list.files(pattern = "csv")
filenames <- filenames[1:3]
for (f in filenames) {
  print(f)
  analyze(f)
}</code></pre>



<div class='out'><pre class='out'><code>[1] "inflammation-01.csv"
</code></pre></div>

<img src="figure/03-loops-R-loop-analyze1.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" /><img src="figure/03-loops-R-loop-analyze2.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" /><img src="figure/03-loops-R-loop-analyze3.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" />

<div class='out'><pre class='out'><code>[1] "inflammation-02.csv"
</code></pre></div>

<img src="figure/03-loops-R-loop-analyze4.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" /><img src="figure/03-loops-R-loop-analyze5.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" /><img src="figure/03-loops-R-loop-analyze6.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" />

<div class='out'><pre class='out'><code>[1] "inflammation-03.csv"
</code></pre></div>

<img src="figure/03-loops-R-loop-analyze7.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" /><img src="figure/03-loops-R-loop-analyze8.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" /><img src="figure/03-loops-R-loop-analyze9.png" title="plot of chunk loop-analyze" alt="plot of chunk loop-analyze" style="display: block; margin: auto;" />

Sure enough, the maxima of these data sets show exactly the same ramp as the first, and their minima show the same staircase structure.

> **Tip:** In this lesson we saw how to use a simple `for` loop to repeat an operation.
As you progress with R, you will learn that there are multiple ways to accomplish this.
Sometimes the choice of one method over another is more a matter of personal style, but other times it can have consequences for the speed of your code.
For instruction on best practices, see this supplementary [lesson](03-supp-loops-in-depth.html) that demonstrates how to properly repeat operations in R.

#### Challenges

1. Write a function called `analyze_all` that takes a filename pattern as its sole argument and runs `analyze` for each file whose name matches the pattern.



<div class="keypoints" markdown="1">
#### Key Points

* Use `for (variable in collection)` to process the elements of a collection one at a time.
* The body of a `for` loop is surrounded by curly braces (`{ }`).
* Use `length(thing)` to determine the length of something that contains other values.
* Use `list.files(pattern = "pattern")` to create a list of files whose names match a pattern.
</div>

#### Next Steps

We have now solved our original problem: we can analyze any number of data files with a single command.
More importantly, we have met two of the most important ideas in programming:

* Use functions to make code easier to re-use and easier to understand.
* Use vectors and data frames to store related values, and loops to repeat operations on them.

We have one more big idea to introduce...
