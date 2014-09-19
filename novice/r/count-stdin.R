lines <- readLines(con = file("stdin"))
count <- length(lines)
cat("lines in standard input: ")
cat(count, sep = "\n")
