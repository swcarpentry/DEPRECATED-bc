# Solutions to exercises for Python Flow Control 
# John Blischak
# jdblischak@gmail.com

"""
Short Exercise
Write an if statement that prints whether x is even or odd.
"""
x = 4
if x % 2 == 0:
    print 'x is even'
else:
    print 'x is odd'

"""
Short Exercise
Using a loop, calculate the factorial of 42 (the product of all integers up to and including 42).
"""
i = 1
factorial = 1
while i <= 42:
    factorial = factorial * i
    i = i + 1
print factorial

"""
Longer Exercise: Converting genotypes

Part 1:
Create a new list which has the converted genotype for each subject ('AA' -> 0, 'AG' -> 1, 'GG' -> 2).
"""
genos = ['AA', 'GG', 'AG', 'AG', 'GG']
genos_new = []
# Use your knowledge of if/else statements and loop structures below.
for i in genos:
    if i == 'AA':
        genos_new.append(0)
    elif i == 'AG':
        genos_new.append(1)
    else:
        genos_new.append(2)

"""
Part 2:
Sometimes there are errors and the genotype cannot be determined. 
Adapt your code from above to deal with this problem (in this example missing data is assigned NA for "Not Available").
"""

genos_w_missing = ['AA', 'NA', 'GG', 'AG', 'AG', 'GG', 'NA']
genos_w_missing_new = []
# The missing data should not be converted to a number, but remain 'NA' in the new list
for i in genos_w_missing:
    if i == 'NA':
        genos_w_missing_new.append(i)
    elif i == 'AA':
        genos_w_missing_new.append(0)
    elif i == 'AG':
        genos_w_missing_new.append(1)
    else:
        genos_w_missing_new.append(2)

"""
Part 3:
The file genos.txt has a column of genotypes. Read in the data and convert the genotypes as above.
Hint: You'll need to use the built-in string method strip to remove the new-line characters 
(See the example of reading in a file above. We will cover string methods in the next section).
"""
# Store the genotypes from genos.txt in this list
genos_from_file = []
handle = open("genos.txt")
for line in handle:
    i = line.strip()
    if i == 'NA':
        genos_from_file.append(i)
    elif i == 'AA':
        genos_from_file.append(0)
    elif i == 'AG':
        genos_from_file.append(1)
    else:
        genos_from_file.append(2)
handle.close()
