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