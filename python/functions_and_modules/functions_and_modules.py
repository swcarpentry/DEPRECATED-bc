# Solutions to exercises for Python Functions and Modules
# John Blischak
# jdblischak@gmail.com

#Short Exercise: Calculate GC content of DNA
# Calculate the fraction of G's and C's in this DNA sequence
seq1 = 'ACGTACGTAGCTAGTAGCTACGTAGCTACGTA'
gc = float(seq1.count('G') + seq1.count('C')) / len(seq1)

"""
Short exercise: Write a function to calculate GC content of DNA
Make a function that calculate the GC content of a given DNA sequence. For the more advanced participants, make your function able to handle sequences of mixed case (see the third test case).
"""

def calculate_gc(x):
    """Calculates the GC content of DNA sequence x.
    x: a string composed only of A's, T's, G's, and C's."""
    x = x.upper()
    return float(x.count('G') + x.count('C')) / (x.count('G') + x.count('C') + x.count('A') + x.count('T'))

"""
Longer exercise: Reading Cochlear implant into Python

Part 1:
Write a function `view_cochlear` that will open the file and print out each line. The only input to the function should be the name of the file as a string. 
"""
def view_cochlear(filename):
    """Reads in data file and prints to console.
    Input: Filename as string.
    """
    x = open(filename)
    for line in x:
        print line.strip()
    x.close()

"""
Part 2:
Adapt your function above to exclude the first line using the flow control techniques we learned in the last lesson. The first line is just `#` (but don't forget to remove the `'\n'`).
"""
def view_cochlear(filename):
    """Reads in data file and prints to console, skipping the first line.
    Input: Filename as string.
    """
    x = open(filename)
    for line in x:
        if line.strip() == '#':
            continue
        else:
            print line.strip()
    x.close()

"""
Part 3:
Adapt your function above to return a dictionary containing the contents of the file. Split each line of the file by a colon followed by a space (': '). The first half of the string should be the key of the dictionary, and the second half should be the value of the dictionary.
"""
def save_cochlear(filename):
    """Reads in data file and saves data.
    Input: Filename as string.
    Output: A dictionary of the data.
    """
    d = {}
    x = open(filename)
    for line in x:
        if line.strip() == '#':
            continue
        else:
            parsed = line.strip().split(': ')
            d[parsed[0]] = parsed[1]
    x.close()
    return d
    
"""
Bonus exercise: Convert DNA to RNA
Write a function that mimics transcription. The input argument is a string that contains the letters A, T, G, and C.
Create a new string following these rules: 
* Convert A to U
* Convert T to A
* Convert G to C
* Convert C to G
Hint: You can iterate through a string using a for loop similary to how you loop through a list.
"""
def transcribe(seq):
    """Transcribes a DNA sequence to RNA.
    Input: string of A's, T's, G's, and C's
    Output: string of RNA basd on input DNA.
    Converts using the following rules:
    A->U, T->A, G->C, C->G
    """
    rna = ''
    for letter in seq:
        if letter == 'A':
            rna = rna + 'U'
        elif letter == 'T':
            rna = rna + 'A'
        elif letter == 'G':
            rna = rna + 'C'
        else:
            rna = rna + 'G'
    return rna
