

# Some best practices for using R and designing programs

1. Start your code with a description of what it is:
	

```r
# This is code to replicate the analyses and figures from my 2014 Science
# paper. Code developed by Sarah Supp, Tracy Teal, and Jon Borelli
```


2. Run all of your import statments (`library` or `require`):


```r
library(ggplot2)
library(reshape)
library(vegan)
```


3. Set your working directory. Avoid changing the working directory once a script is underway. Use `setwd()` first . Do it at the beginning of a R session. Better yet, start R inside a project folder.

4. Use `#` or `#-` to set off sections of your code so you can easily scroll through it and find things.

5. If you have only one or a few functions, put them at the top of your code, so they are among the first things run. If you written many functions, put them all in their own .R file, and `source` them. Source will run all of these functions so that you can use them as you need them.


```r
source("my_genius_fxns.R")
```


6. Use consistent style within your code. 

7. Keep your code modular. If a single function or loop gets too long, consider breaking it into smaller pieces.

8. Don't repeat yourself. Automate! If you are repeating the same piece of code on multiple objects or files, use a loop or a function to do the same thing. The more you repeat yourself, the more likely you are to make a mistake.

9. Manage all of your source files for a project in the same directory. Then use relative paths as necessary. For example, use


```r
dat <- read.csv(file = "/files/dataset-2013-01.csv", header = TRUE)
```


rather than:


```r
dat <- read.csv(file = "/Users/Karthik/Documents/sannic-project/files/dataset-2013-01.csv", 
    header = TRUE)
```


10. Don't save a session history (the default option in R, when it asks if you want an `RData` file). Instead, start in a clean environment so that older objects don't contaminate your current environment. This can lead to unexpected results, especially if the code were to be run on someone else's machine.

11. Where possible keep track of `sessionInfo()` somewhere in your project folder. Session information is invaluable since it captures all of the packages used in the current project. If a newer version of a project changes the way a function behaves, you can always go back and reinstall the version that worked (Note: At least on CRAN all older versions of packages are permanently archived).

12. Collaborate. Grab a buddy and practice "code review". We do it for methods and papers, why not code? Our code is a major scientific product and the result of a lot of hard work!

13. Develop your code using version control and frequent updates!

### Challenges

1. What other suggestions do you have?
2. How could we restructure the code we worked on today, to make it easier to read? Discsuss with your neighbor.
3. Make two new R scripts called inflammation.R and inflammation_fxns.R 
4. Copy and paste the code so that inflammation.R "does stuff" and inflammation_fxns.R holds all of your functions. __Hint__: you will need to add `source` code to one of the files.

