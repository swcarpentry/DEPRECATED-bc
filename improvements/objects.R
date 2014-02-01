# Objects in R

## type or mode: how the object is stored
mode()

## class or data structure: what info it contains and how to use the object
class()

## VECTOR
(v1 <- c(1, 3, 5))
(v2 <- c(1, 3, "a")) # coercion
mode(v2)
as.character(v1)
1:10
seq(from = 5, to = 25, by = 5)

### FACTOR
(breed <- factor(c("Holstein", "Brown Swiss", "Holstein", "Ayrshire",
                   "Canadian")))
table(breed)

## LIST
(l <- list("a", "b", "c"))
(cow <- list(breed = "Holstein", age = 3, last.prod = c(25, 35, 32)))
cow$breed
cow[[1]]
(h <- hist(islands))
str(h)

## MATRIX
(m <- rbind(c(1, 4), c(2, 2)))
(m <- matrix(data = 1:12,
             nrow = 4, ncol = 3,
             dimnames = list(c("cow1", "cow2", "cow3", "cow4"),
                 c("milk", "fat", "prot"))))

## ARRAY
(a <- array(data = 1:24, dim = c(3, 4, 2)))

## DATA FRAME
(df <- data.frame(cow = c("Moo-Moo", "Daisy", "Elsie"),
                  prod = c(35, 40, 28),
                  pregnant = c(TRUE, FALSE, TRUE)))

