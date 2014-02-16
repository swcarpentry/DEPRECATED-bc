
# Exercise 02

1. Take the `iris` dataset. Using `Species` as the indicator variable, convert the wide dataset to long such that it appears as below:

```
  Species     variable value
1  setosa Sepal.Length   5.1
2  setosa Sepal.Length   4.9
3  setosa Sepal.Length   4.7
4  setosa Sepal.Length   4.6
5  setosa Sepal.Length   5.0
6  setosa Sepal.Length   5.4
```

2.  Save this dataset to disk as `"species1.csv"`. Duplicate the previous function call and save it again as `"species2.csv"``

Write code that does the following:

Create a list of length two containing "species1.csv" and "species2.csv". 
Now write a `lapply` function to read both csv files using write.csv into a list.

__Advanced version__
* Reads all file names in your working directory with a `.csv` extension.  
* Treating that as a list, run the `read.csv` function on each item so you get a list of size 2 that contains all the data.

*Note: If this is a hard problem. Work with your neighbor*

3. With the list from the previous result, write a `lappy or llply` function call that will retrieve row number 1 from each dataset and return that back to a list.

4. Read the `gapminderDataFiveYear.txt` dataset in the data folder into an object. Split the data by continent and country, find the year with the highest life expectancy for each combination and return those results back into a data.frame. 

Hint: You will write a anonymous function (an unnamed function) that will exist only inside one of plyr's function call. 

e.g. 
```
ddply(iris, .(Species), function(x) {
    x[1, ]
})
# This trivial function returns the first row of each species. Note that I did not name this function.
```

You will have to:
    split by continent and country
    Inside your anonymous function find the year with the highest life expectancy
    return that record back.

Advanced version of this question: Also return the year with highest population and highest per capita gdp. So you'll return 3 rows per country. 


5. Read the `mammals.csv` file into a data.frame. Use `ddply` to split the dataset by limb morphology, then write each file to a separate text file named by limb.

Hints: Use `unique` inside your anoynymous function to get a unique name.
Use `paste` to create a filename. e.g. `paste(unique(x$Limb_morphology), ".csv", sep="")`

Pass this variable to the `file` argument in `write.csv`



