# Literate programming in R using `knitr`

It's easy to generate reports dynamically in `R`. Literate programming is a paradigm that originally came from Donald Knuth (Knuth, 1984).

**Basic idea:** Write **data** + **software** + **documentation** (or in this case manuscripts, reports) together.

Analysis code can be divided into text and code "chunks". Doing so allows us to extract the code for machine readable documents (technically referred to as a `tangle`) or produce a human-readable document (also called `weave`).

Literate programming involves three main steps:  

1. Separate the narrative from the code
2. Execute source code and return the results.
3. Combine the results from the source code with the original narratives to produce a final document.

## Why this is important?
Results from scientific research have to be easy to reproduce so others can verify results making them trustworthy. Otherwise we risk producing one off results that no one outside the original research group can reproduce. In this lesson we will learn reproducible research, which is one of the by products of dynamics report generation. However, this process alone will not always guarantee reproducibility. 


## Installing `knitr`

```coffee
# Installing knitr is quite easy. 
install.packages("knitr", dependencies = TRUE)
install.packages("pander") # Pander is a useful package for formatting tables in markdown.
```

Knitr supports a variety of documentation formats including `markdown`, `html` and `LaTeX`. It also allows for easy export to `PDF` and `HTML`.

## What is markdown?

Markdown is an incredibly simple semantic file format, not too dissimilar from .doc, .rtf, or .txt. Markdown makes it easy for even those without significant knowledge of markup languages like html or LaTex to write any sort of text (including with links, lists, bullets, etc.) and have it parsed into a variety of formats. To learn more about the basics of markdown, peruse this [short tutorial on the format](markdown.md).

* [Markdown cheatsheet](https://github.com/adam-p/markdown-here/wiki/Markdown-Cheatsheet)
* [Original markdown reference](http://daringfireball.net/projects/markdown/basics)


## When to implement reproducibility via literate programming?

* You can do this anytime but it is always best to do this at the beginning of a project. Much like writing unit tests, waiting till after a project is completed is a bad idea. 
* Best used alongside version control like Git. 
* Use software like R where instructions are scripted.
* Never save the output (only raw dataset with pre-processing code. Even cleaned datasets can be discarded but it might help to save them temporarily during intermediate steps).
* Finally, store the data in a non-proporietary format (e.g. `csv` over `xls`). This will ensure that your data will be readable long into the future. 


---

**Good and Bad Practices**

Related: See [best practices](../R-basics/best-practices.Rmd) in the [R-basics folder](../R-basics/).


## Creating a basic knitr document

In RStudio, choose new R Markdown file (easiest way)
or you can create a new text file and save it with extension `.Rmd`.

A basic code chunk looks like this:

<pre><code>```{r}
# some R code
```
</code></pre>
---

You can knit this document using the knit button or do it programmatically using the `knit()` function.

```
library(knitr)
knit("file.Rmd")
```

**What just happened?**

knitr read the Rmd file, then located and ran all the code chunks identified by the backticks, and replaced it with the output of R function calls. If figures are generated from any such calls, they will be included in markdown syntax.  


## Chunk labels

You can also name your code chunks. This allows you to keep all the code in a separate script and just refer to code chunks using meaningful names (e.g. data-processing, analysis, model-fitting, visualization, figures, tables)


```
{r, chunk_name}
```


__Some rules on naming chunks__
* Chunk labels are supposed to be unique id’s in a document.  
* knitr will throw an error if two chunks have the same name.  
* If no chunk names are given, knitr will simply increment from chunk 1, 2,3 etc.

In addition to naming chunks within the curly braces, you can also add a bunch of other options on how that particular code chunk should behave. 

**Other options you can add to the tag**

| Option | Description |
| ------ | ------------ |
| **echo** =   TRUE or FALSE |  to show or hide code.  |
| **eval** =   TRUE or FALSE | to run or skip the code.  |
| **warning** =   TRUE or FALSE | to show or hide function warnings.  |
| **message** =  TRUE or FALSE | to show or hide function R messages.  |
| **results** = "hide"  | will hide results. They will still be executed |
| **fig.height** = | Height of figure  |
| **fig.width** =  | width of figure  |

Once your output markdown files (`.md`) files are generated, you should never edit them because they are automatically generated. Next time you knit the original `.Rmd` files, all the changes in the `.md` file will get wiped out. 

**Write sentences in text with inline output**

```
Include some text `r mean(1:5)`. 
```

**Summarizing output from models.**

<pre><code>```{r fit_model}
library(datasets)
data(airquality)
fit = lm(Ozone ~ Wind + Temp + Solar.R, data = airquality)
```

## Including formatted tables in markdown

```{r showtable, results="asis", echo = FALSE, message = FALSE, eval = FALSE, warning = FALSE}
library(pandoc)
pander(fit)
```</code></pre>

## Global options

Global options are shared across all the following chunks after the location in which the options are set, and local options in the chunk header can override global options.

<pre><code>```{r setoptions, eval = FALSE, echo = FALSE}
options(width = 60, show.signif.stars = FALSE)
opts_chunk$set(echo = FALSE, 
            results = "asis", 
            warning = FALSE, 
            message = FALSE, 
            fig.width = 5, 
            fig.height = 4, 
            tidy = TRUE, 
            fig.align = 'center')
```</code></pre>


## Other Options
**Dealing with long running process**

By adding `cache = TRUE` to a code block definition. After the first run, results will be cached. We'll discuss better ways to acheive the same thing using `Make` in the next section.


## Quick reporting

Generating reports in knitr doesn't always have to involve a laborious `.Rmd ` file where scripts need to be broken down into smaller chunks. Sometimes a user might need a simple report generated very quickly from an existing script. The function `stitch()` in knitr makes it possible to generate nicely formatted reports from R scripts. 

knitr provides a template of the source document with default settings which  allows the user to simply pass any R script into this template (consider this one giant code chunk). knitr will  then compile the template to a report. The package currently has build in support for a range of templates from LateX, html, and markdown. To stitch a report:

```
library(knitr) 
stitch("your-script.R")
```
## Additional chunk options

Chunks are extremely flexible and more options (beyond the ones listed in the table above) can be included in the header. These look exactly like the kinds of arguments that one might pass to standard R functions. In the example below, the chunk will only be executed if the condition (in this case x less than 5) is satisfied. 

```
{r chunk_name, eval = if (x < 5) TRUE else FALSE}
```
This allows your document to be dynamic allowing certain chunks to be executed only when specific conditions are met.



## Error handling

By default `knitr` will not stop execution if it encounters an error. It will continue through to the end of the document and include any errors that arise within chunks. The reason for this behavior is that knitr treats the code as if it were fed directly into the R console. Any errors get printed to the screen and the remaining commands are executed. 

To stop knitr as soon as it encounters an error, one can set an option explicitly:

```coffee
opts_knit$set(stop_on_error = 2L)
```

__Possible options for error handling__

| Option | What it does | 
| ------ | ------------ | 
| * **0L** | do not stop on errors, continue on as if code was pasted into R console  |
| * **1L** | when an error occurs, return the results up to this point, ignore the rest of the code within that particular chunk without reporting any further errors. |
| * **2L** | Completely stop upon encountering the first error. |


## Working with graphics in `knitr`

If you use `ggplot2` from the data visualization section, you can have that easily parsed into your document.

```
library(ggplot2)
p <- qplot(carat, price, data = diamonds) + geom_hex()
p 
# no need to explicitly print(p)
```

---

## Code Externalization

Sometimes it can be rather tedious to include dozens of lines of code in the same file as the narrative. In such cases, one can improve readability by externalizing the code into a separate script and simply calling the chunks at the appropriate locations in a document. The code will get read in and executed at those points. There are two steps to making this happen.

First, read the script using `read_chunk()` at the top of any `.Rmd` file
```
read_chunk("source_code.R")
```

This is usually done in an early chunk such as the first chunk of a document, and we can use the chunk `data-processing` later in the source document:


```
{r, data-processing}
```

Then simply call any chunk as needed simply by using its label. You do not have to include any code between the backticks. 

```
## @knitr data-processing
```

Be sure to leave a blank line between chunks.

--- 

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

You can add more options to the basic Pandoc call. To see a full list of options

```
pandoc --help
```

A commonly used option is to add margins using the `-V` argument (in this case 1 inch):

```
pandoc -V geometry:margin=1in test.md -o test.pdf
```

You can use `Make` to automate much of this process. For example, by setting up a series of dependencies, you can have a new document knitted if the underlying data or code changes but not otherwise. 
