count <- 0
lines <- readLines(con = file("stdin"))
for (line in lines) {
  count <- count + 1
}

cat("lines in standard input: ")
cat(count)