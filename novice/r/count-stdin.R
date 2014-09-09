count <- 0
lines <- readLines(con = file("stdin"))
for (line in lines) {
  count <- count + 1
}

print("lines in standard input:")
print(count)