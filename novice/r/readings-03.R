main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  for (filename in args) {
    dat <- read.csv(file = filename, header = FALSE)
    mean_per_patient <- apply(dat, 1, mean)
    print(mean_per_patient)
  }
}

main()