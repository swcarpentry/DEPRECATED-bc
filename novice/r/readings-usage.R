main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  action <- args[1]
  filenames <- args[-1]
  if (!(action %in% c("--min", "--mean", "--max"))) {
    usage()
  } else if (length(filenames) == 0) {
    process(file("stdin"), action)
  } else {  
    for (f in filenames) {
      process(f, action)
    }
  }
}

process <- function(filename, action) {
  dat <- read.csv(file = filename, header = FALSE)
  
  if (action == "--min") {
    values <- apply(dat, 1, min)
  } else if (action == "--mean") {
    values <- apply(dat, 1, mean)
  } else if (action == "--max") {
    values <- apply(dat, 1, max)
  }
  cat(values, sep = "\n")
}

usage <- function() {
  cat("usage: Rscript readings-usage.R [--min, --mean, --max] filenames", sep = "\n")
}

main()