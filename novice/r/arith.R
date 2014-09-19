main <- function() {
  # Performs addition or subtraction from the command line.
  #
  # Takes three arguments:
  # The first and third are the numbers.
  # The second is either + for addition or - for subtraction.
  #
  # Ex. usage:
  #   Rscript arith.R 1 + 2
  #   Rscript arith.R 3 - 4
  #
  args <- commandArgs(trailingOnly = TRUE)
  num1 <- as.numeric(args[1])
  operation <- args[2]
  num2 <- as.numeric(args[3])
  if (operation == "+") {
    answer <- num1 + num2
    cat(answer)
  } else if (operation == "-") {
    answer <- num1 - num2
    cat(answer)
  } else {
    stop("Invalid input. Use + for addition or - for subtraction.")
  }
}

main()