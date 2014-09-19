main <- function() {
  # Checks that all csv files have the same number of rows and columns.
  #
  # Takes multiple arguments: the names of the files to be checked.
  #
  # Ex. usage:
  #   Rscript check.R inflammation-*
  #
  args <- commandArgs(trailingOnly = TRUE)
  first_file <- read.csv(args[1], header = FALSE)
  first_dim <- dim(first_file)
#   num_rows <- dim(args[1])[1]  # nrow(args[1])
#   num_cols <- dim(args[1])[2]  # ncol(args[1])
  for (filename in args[-1]) {
    new_file <- read.csv(filename, header = FALSE)
    new_dim <- dim(new_file)
    if (new_dim[1] != first_dim[1] | new_dim[2] != first_dim[2]) {
      cat("Not all the data files have the same dimensions.")
    }
  }
}

main()