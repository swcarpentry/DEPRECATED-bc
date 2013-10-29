# Literate programming in R

It's easy to generate reports dynamically in R. The paradigm come from Donald Knuth (Knuth, 1984).

Basic idea: Write software + documentation (or in this case manuscripts, reports) together.

This allows us to extract the code (technically called `tangle`) or produce a document (called `weave`).

Literate programming involves with three main steps:
1. parse the source document and separate code from narratives
2. execute source code and return results
3. mix results from the source code with the original narratives


## What is markdown?

An incredibly simple semantic file format, not too dissimilar from .doc, .rtf or .txt. Markdown makes it easy for even those without a web-publishing background to write prose (including with links, lists, bullets, etc.) and have it displayed like a website. 

* [Markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
* [Original markdown reference](http://daringfireball.net/projects/markdown/basics)


## Installing Knitr


```coffee
install.packages("knitr", dependencies = TRUE)
```



----
# OLD NOTES
# How to use knitr

* The basic idea of dynamic report generation stems from literate programming, a programming paradigm conceived by Donald Knuth (Knuth, 1984). The original idea was mainly for writing software: mix the source code and documentation together; we can either extract the source code out (called tangle) or execute the code to get the compiled results (called weave).
* Technically, literate programming involves with three steps:
1. parse the source document and separate code from narratives
2. execute source code and return results
3. mix results from the source code with the original narratives

Results from scientific research have to be reproducible to be trustwor- thy. We do not want a finding to be merely due to an isolated occur- rence, e.g. only one specific laboratorian can produce the results on one specific day, and nobody else can produce the same results under the same conditions.

Reproducible research (RR) is one possible by-product of dynamic report generation, but the latter does not absolutely guarantee RR. Be- cause there is usually no human intervention when we generate a re- port dynamically, it is likely to be reproducible since it is relatively easy to prepare the same software and hardware environment, which is ev- erything we need to reproduce the results. 

**Good and Bad Practices**

The key to keep in mind for RR is that other people should be able to reproduce our results, therefore we should try out best to make out computation portable. We discuss some good practices for RR below and explain why it can be bad not to follow them.

• manage all source files under the same directory and use relative paths whenever possible: absolute paths can break reproducibility, e.g. a data file like C:/Users/someone/foo.csv or /home/someone/foo.csv may only exist in one computer, and other people may not be able to read it since the absolute path has changed in their hard disk; if we keep everything under the same directory, we can use read a data file like read.csv(’foo.csv’) (if it is under the current working direc- tory) or read.csv(’../data/foo.csv’) (go one upper level and find the file under the data/ directory); when we disseminate the results, we can make an archive of the whole directory (e.g. as a zip package)

• do not change the working directory after the computing has started: setwd() is the function in R to set the working directory, and it is not uncommon to see setwd(’path/to/some/dir’) in user’s code, which is bad because it is not only an absolute path, but also has a global effect to the rest of the source document – we have to keep in mind that all relative paths may need adjustments since the root directory has changed, and the software may write the output in an unexpected place (e.g. the figures are expected to be generated to ./figures/ directory, but are actually written to ./data/figures/ instead if we setwd(’./data/’)); if we have to set the working directory at all, do it in the very beginning of an R session; most of the editors to be introduced in Chapter 4 follow this rule, and the working directory is set to the directory of the source document before knitr is called

• compile the documents in a clean R session: existing R objects in the current R session may “contaminate” the results in the output; it is fine that we write a report by accumulating code chunks one by one and running them interactively to check the results, but in the end we should compile a report in the batch mode with a new R session so all the results are freshly generated from the code
￼￼￼￼￼
• avoid environment variables for data analysis: while environment variables are often heavily used in programming for configuration purposes, it is ill-advised to use them in data analysis because they require additional instructions for users to set up, and humans can simply forget to do this; if we have any options to set up, do it inside the source document

• attach: sessionInfo()and instructions on how to compile this document: the session information makes a reader aware of the software environment such as the version of R, the operating system and add-on packages used; sometimes it is not as simple as calling one single function to compile a document, and we have to make it clear how to compile it if additional steps are required, but it is better to provide the instructions in the form of a computer script, e.g. a Shell script, a Makefile or a batch file.


## How to install the package:

```coffee
install.packages("knitr", dependencies = TRUE)
```


## Options

Some options for figures

e.g.
```
{r model, fig.width=4, fig.height=3, fig.align='center'}
```


The ideal output from Markdown is an HTML web page, as shown in Figure 3.4 (in Mozilla Firefox). Similarly we can see the syntax for R code in a Markdown document: ```{r} opens a code chunk, ``` terminates a chunk and inline R code can be put inside`r `,where`is a backtick.

## Quick reporting

Quick Reporting
If a user only has basic knowledge about R but knows nothing about knitr, or one does not want to write anything else other than an R script, it is also possible to generate a quick report from this R script using the stitch() function.

The basic idea of stitch() is that knitr provides a template of the source document with some default settings, so that the user only needs to feed this template with an R script (as one code chunk), then knitr will compile the template to a report. Currently it has built-in templates for LATEX, HTML and Markdown. The usage is like this:

```coffee
library(knitr) stitch("your-script.R")
```
## Chunk options

As mentioned in Chapter 3, we can write chunk options in the chunk header. The syntax for chunk options is almost exactly the same as the syntax for function arguments in R. They are of the from
option = value
There is nothing to remember about this syntax due to the con- sistency with the syntax of R: as long the option values are valid R code, they are valid to knitr. Besides constant values like echo = TRUE (a logical value) or out.width = ’\\linewidth’ (character string) or fig.height = 5 (a number), we can write arbitrary valid R code for chunk options, which makes a source document programmable. Here is a trivial example:

```
<<foo, eval=if (bar < 5) TRUE else FALSE>>=
```

Suppose bar is a numeric variable created in the source document before this chunk.We can pass an expression if (bar < 5) TRUE else FALSE to the option eval, which makes the option eval depends on the value of bar, and the consequence is we evaluate this chunk based on the value of bar (if it is greater than 5, the chunk will not be evaluated), i.e. we are able to selectively evaluate certain chunks. 


## Chunk Label
The only possible exception is the chunk label, which does not have to follow the syntax rule. In other words, it can be invalid R code. This is for both historical reasons (Sweave convention) and laziness (avoid typing quotes). Strictly speaking, the chunk label, as a part of chunk options, should take a character value, hence it should be quoted, but in most cases, knitr can take care of the unquoted labels and quote them internally, even if the “objects” used in the label expression do not exist. Below are all valid ways to write chunk labels:

```
<<echo=FALSE, label="foo-bar">>=
```

Chunk labels are supposed to be unique id’s in a document, and they are mainly used to generate external files such as images (Chapter 7) and cache files (Chapter 8). If two non-empty chunks have the same label, knitr will stop and emit an error message, because there is poten- tial danger that the files generated from one chunk may override the other chunk. If we leave a chunk label empty, knitr will automatically generate a label of the form unnamed-chunk-i, where i is an incremen- tal chunk number from 1, 2, 3,·

## Global Options
Chunk options control every aspect of a code chunk, as we will see in more details through Chapter 6 to 11. If there are certain options which are used commonly for most chunks, we can set them as global chunk options using the object opts chunk. Global options are shared across all the following chunks after the location in which the options are set, and local options in the chunk header can override global options. For example, we set the option echo to FALSE globally:

```
opts_chunk$set(echo = FALSE, message = FALSE)
```

# Including inline code in markdown

```
`r x`
```

Inline math: $\alpha + \beta$. Display style:
$$f(x) = xˆ2 + 1$$

the inline hook is not associated with a code chunk; it defines how to format the output from inline R code, for example, we may want to round all the numbers from inline output to 2 digits and we can define the inline hook as:
in fact knitr takes care of rounding in the default inline hook (Section 6.1), so we do not really have to reset this hook;

```
knit_hooks$set(inline = function(x) { if (is.numeric(x))
x <- round(x, 2)
as.character(x) # convert x to character and return
})
```


Sometimes we do not want to mix R code with normal texts, but write texts in comments, so that the whole document is a valid R script. The function spin() in knitr can deal with such R scripts if the comments are written using the roxygen syntax. The basic idea of spin() is still lit- erate programming: when we compile this R script, #' will be removed so that normal texts are “restored”, and R code will be evaluated. Any- thing that is not behind a roxygen comment is treated as a code chunk. To write chunk options, we can use another type of special comments #+ or #- followed by chunk options. Below is a simple example:
#' Introduce the method here; then write R code: 1+1
x <- rnorm(10)
#' It is also possible to write chunk options, e.g.
#'
#+ test-label, fig.height=4
plot(x)
#' The document is done now.
We can save this script to a file called test.R, and compile it to a report:
The spin() function has a format argument which specifies the out- put document format (default to R Markdown). For example, if format = ’Rnw’, the R code will be first inserted between <<>>= and @, and then compiled to generate LATEX output.
This looks similar to the stitch() function in Section 3.3, which also creates a report based on an R script, but spin() makes it possible to write text chunks and stitch() can only use a predefined template, so there is less freedom.


Another R option digits controls how many digits a number should be rounded to; knitr uses 4 by default, and R’s default is 7 which of- ten makes a number unnecessarily “precise”. For example, a number 123456789 will become 1.2346 × 108 in the final output. We can change the defaults in the first chunk of a document like:

```
## numbers >= 10ˆ5 will be denoted in scientific 
## notation, and rounded to 2 digits 
options(scipen = 1, digits = 2)
```


The chunk option eval (TRUE or FALSE) decides whether a code chunk should be evaluated. When a chunk is not evaluated, there will be no results returned except the original source code. This option can also take a numeric vector to specify which expressions (when chunk op- tion tidy = TRUE) or lines (when tidy = FALSE) to be evaluated; in this case, the code that is set not to be evaluated will be commented out. 


Note these two options are not specific to knitr; they are global op- tions in R. If we are not satisfied with the default inline output, we can rewrite the inline hook as introduced in Section 5.3. 


Syntax highlighting comes by default in knitr (chunk option highlight = TRUE) since it enhances the readability of the source code – character strings, comments and function names, etc, are in different colors. This option only works for LATEX and HTML output, and it is not necessary for Markdown because there are other libraries which can highlight code in web pages, e.g. RStudio uses a JavaScript library highlight.js to do syntax highlighting for Markdown output.

## This is good.

```
{r, tidy = TRUE}
```

The function tidy.source() in the formatR package (Xie, 2012a) is used to reformat R code (option tidy = TRUE), e.g. it can add spaces and in- dentation, break long lines into shorter ones and automatically replace the assignment operator = to <-; 


## Show/Hide Output
We can show or hide different parts of the text output including the source code, normal text output, warnings, messages, errors and the whole chunk. Below are the corresponding chunk options with default values in the braces:
echo (TRUE)whethertoshowthesourcecode;itcanalsotakeanumeric vector like the eval option to select which expressions to show in the output, e.g. echo = 1:3 selects the first 3 expressions, and echo = -5 means do not show the 5th expression
results (’markup’)howtowrapupthenormaltextoutputwhichwould have been printed in the R console if we run the code in R; the default value means to mark up the results in special environments such as LATEX environments or HTML div tags; two other possible values are
’asis’ write the raw output from R to the output document without any markups, e.g. the source code cat(’<em>emphasize</em>’) can produce an italic text in HTML when results = ’asis’; this is very useful when we use R to produce raw elements for the output, e.g. tables
’hide’ this option value hides the normal text output
warning/error/message (TRUE) whether to show warnings, errors and messages; usually these three types of messages are produced by warn- ing(), stop() and message() in R
split (FALSE) whether to redirect the chunk output to a separate file (the filename is determined by the chunk label); for LATEX, \input{} will be used if split = TRUE to input the chunk output from the file; for HTML, the <iframe> tag will be used; other output formats will ignore this option
include (TRUE) whether to include the chunk output in the document; when it is FALSE, the whole chunk will be absent in the output, but the code chunk will still be evaluated unless eval = FALSE


As we have introduced in Section 5.1, we can use opts chunk to set global chunk options. For instance, we want to suppress all warnings and messages in the whole document, then we can do this in the first chunk of the document:

```
opts_chunk$set(warning = FALSE, message = FALSE)
```

It may be very surprising to knitr users that knitr does not stop on errors! As we can see from the previous example, 1 + ’a’ should have stopped R because that is not a valid addition operation in R (a number + a string). The default behavior of knitr is to act as if the code were pasted into an R console: if you paste 1 + ’a’ to the R console, you will see an error message, but that does not halt R – you can continue to type or paste more code. To completely stop knitr when errors occur, set this option in advance:

```
opts_knit$set(stop_on_error = 2L)
```


The meaning of the integer code for stop on error is as follows
(from the evaluate package; see the documentation for evaluate() there):
0L do not stop on errors, just like the code was pasted into R console
1L when an error occurs, return the results up to this point and ignore the rest of code in the chunk but do not throw the error either
2L a full stop on errors

## Themes

See knitr themes

```
head(knit_theme$get(), 20)
```

```
knit_theme$set("autumn")
```

knitr themes: http://animation.r-forge.r-project.org/knitr/

Note syntax highlighting themes only work for LATEX and HTML output. For Markdown, the highlight.js library also allows customiza- tion but that is beyond the scope of R and knitr. See http://bit.ly/ knitr-themes for a preview of all these themes.

## Working with graphics in knitr

```coffee
￼library(ggplot2)
p <- qplot(carat, price, data = diamonds) + geom_hex() p # no need to print(p)
```


# Plot Size in Output
The fig.width and fig.height options specify the size of plots in the graphical device, and the real size in the output document can be dif- ferent (specified by out.width and out.height). When there are mul- tiple plots per code chunk, it is possible to arrange multiple plots side by side.


To solve the second problem, we need to let knitr know changes in external files. One natural indicator is the modification time of files, which can be obtained by the function file.info(). Suppose the data file is named iris.csv, and we can put its modification time in a chunk option iris time, e.g.

```
<<cache-rversion, cache=TRUE, version=R.version.string>>=
# code which may be affected by R version
R.version.string
## [1] "R version 2.15.2 (2012-10-26)"
@
￼￼<<itime, cache=TRUE, iris_time=file.info('iris.csv')$mtime>>
# data will be re-read if iris.csv becomes newer
iris <- read.csv("iris.csv") @
```

## More on the setup chunk

```
<setup, cache=FALSE, include=FALSE>>=
# set up some global options for the document
options(width = 60, show.signif.stars = FALSE)
# also set up global chunk options
library(knitr)
opts_chunk$set(fig.width = 5, fig.height = 4, tidy = FALSE) @
```

### Manual dependency
<<chunkB, dependson='chunkA'>>=


## Code Externalization
It can be more convenient to write R code chunks in a separate R script, rather than mixing them into a source document; 

Labeled Chunks
The setting is like this: the R script also uses chunk labels (marked in theform## @knitr chunk-label);ifthecodechunkinthesourcedoc- ument is empty, knitr will match its label with the label in the R script to input external R code.
For example, suppose this is a code chunk labelled as Q1 in an R script named shared.R which is under the same directory as the source document:
￼
```
## @knitr Q1
gcd <- function(m, n) {
while ((r <- m%%n) != 0) {
m <- n n <- r
}n }
```

In the source document, we can first read the script using the func- tion read chunk():

```
read_chunk("shared.R")
```


This is usually done in an early chunk such as the first chunk of a document, and we can use the chunk Q1 later in the source document:


```￼
￼<<Q1>>
= @
```


Hooks are an important component to extend knitr. A hook is a user- defined R function to fulfill tasks beyond the default capability of knitr. There are two types of hooks: chunk hooks and output hooks. We have already introduced some built-in output hooks in Section 5.3, and how to customize both the chunk and inline R output. In this chapter we focus on chunk hooks.

```
￼names(knit_hooks$get(default = T))
```

## Tricks and Solutions
In this chapter we show some tricks which can be useful for writing and compiling reports more easily and quickly, and also solutions to frequently asked questions.


Option Aliases
We may feel some options are very frequently used but the names are too long to type. In this case we can set up aliases for chunk options using the function set alias() in the beginning of a document, e.g.

```
set_alias(w = "fig.width", h = "fig.height")
```

Then we will be able use w and h for the figure width and height
￼￼￼￼
respectively, e.g.
The chunk above is equivalent to:
```
￼<<fig-size, w=5, h=3>>= 
plot(1:10)
@
```

```
opts_template$set(
fig.large = list(fig.width = 7, fig.height = 5), 
fig.small = list(fig.width = 3.5, fig.height = 3)
)
```

```
￼<<fig-ex, opts.label='fig.large'>>= 
plot(1:10)
@
```


## Pandoc
Pandoc (http://johnmacfarlane.net/pandoc) is a universal document converter. In particular, Pandoc can convert Markdown to many other document formats, including LATEX, HTML, Rich Text Format (*.rtf), E- Book (*.epub), Microsoft Word (*.docx) and OpenDocument Text (*.odt), etc.
Pandoc is a command line tool. Linux users should be fine with it; for Windows users, the command window can be accessed via the Start menu, then Run cmd. Once we have opened a command window (or terminal), we can type commands like this to convert a Markdown file, say, test.md, to other formats:

```
pandoc test.md -o test.html
pandoc test.md -o test.odt
pandoc test.md -o test.rtf
pandoc test.md -o test.docx
pandoc test.md -o test.pdf
```
