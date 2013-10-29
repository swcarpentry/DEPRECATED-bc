
# Input output operations

## Inputting data

```coffee
x <- scan("data_file.txt")
# add a separator
x <- scan("data/messy_data.txt", what=" ", sep = "\n")
# or read data from the console
x <- scan()
# keep entering values and hit an empty return key to end
```
Reading single lines (e.g. user input)

```coffee
variable <- readline()
# or provide more information
variable <- readline("Enter number of simulations: ")
```


## Reading files  
Most plain text files can be read with `read.table` or variants thereof (such as `read.csv`).

```coffee
df <- read.table("data.dat", header = TRUE)
```

or using `readLines`

```coffee
dt <- readLines("data/messy-data.txt")
```

## Files from the web

```
url <- "https://raw.github.com/karthikram/ggplot-lecture/master/climate.csv"
my_data <- read.csv(url, header = TRUE)
```

## Local file operations

One can list files from any local source as well.

```coffee
list.files()
file.info()
dir()
file.exists()
getwd()
setwd()
```


---



## Writing files

Saving files is easy in R. We have loaded the `iris` dataset into our memory. Can you save this back to a `csv` file to disk with the name `tgac_iris.csv`?

What commands did you use?


# Short term storage

```coffee
saveRDS(iris, file = "tgac_iris.rds")
iris_data <- readRDS("tgac_iris.rds")
```
This is great for short term storage. All factors and other modfications to the dataset will be preserved. However, only R can read these data back and not the best option if you want to keep the file stored in the easiest format.

# Long-term storage

```coffee
write.csv(iris, file = "tgac_iris.csv", row.names = FALSE)
```

![](saving_files.png)

# Easy to store compressed files to save space:

```coffee
write.csv(diamonds, file = bzfile("diamonds.csv.bz2"),
  row.names = FALSE)
```

# Reading is even easier:

```coffee
diamonds5 <- read.csv("diamonds.csv.bz2")
```

Files stored with `saveRDS()` are automatically compressed.
