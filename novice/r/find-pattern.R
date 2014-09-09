main <- function() {
  # Finds all files in the current directory that contain a given pattern.
  #
  # Takes one argument: the pattern to be searched.
  #
  # Ex. usage:
  #   Rscript find-pattern.R csv
  #
  args <- commandArgs(trailingOnly = TRUE)
  pattern <- args[1]
  files <- list.files(pattern = pattern)
  cat(files, sep = "\n")
}

main()