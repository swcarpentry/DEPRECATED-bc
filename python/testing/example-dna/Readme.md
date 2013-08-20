**Materials originally by John Blischak**

## A TDD Example

To illustrate TDD, let's return to the function you wrote yesterday,
`calculate_gc`. We'll start from scratch and develop the function
by meeting test specifications. 

The beginning of the function is contained in the file `calculate_gc.py`
in this directory. It currently takes one argument as input, but does
nothing.

```python
def calculate_gc(x):
    '''
    Calculates the GC content of DNA sequence x.
    '''
    pass
```

The tests that we must pass are contained in the file
`test_calculate_gc.py`. We can run the tests using `nosetests`.

    nosetests -v test_calculate_gc.py

As expected, we fail all the tests! What is the bare minimum 
functionality we must add to pass the first test below?

```python
def test_only_G_and_C():
    '''
    Sequence of only G's and C's has fraction 1.0
    '''
    fixture = 'GGCGCCGGC'
    result = calculate_gc(fixture)
    assert_equal(result, 1.0)
```

And the second test?

```python
def test_half():
    '''
    Sequence with half G and C has fraction 0.5
    '''
    fixture = 'ATGC'
    result = calculate_gc(fixture)
    assert_equal(result, 0.5)
```

Test number three?

```python
def test_lower_case():
    '''
    Sequence with lower case letters
    '''
    fixture = 'atgc'
    result = calculate_gc(fixture)
    assert_equal(result, 0.5)
```

Test number four?

```python
def test_not_DNA():
    '''
    Raise TypeError if not DNA
    '''
    fixture = 'qwerty'
    assert_raises(TypeError, calculate_gc, fixture)
```

Through this cycle of writing tests and modifying the function to pass 
the tests, we have developed a function that behaves exactly as we 
expect and nothing more. And the tests not only serve as documentation 
of what the function does, but can also be easily ran again if we made 
further modifications (regression tests). What would be the next test 
you would write for our function?

## Exercise: Test function that transcribes DNA to RNA

In the lesson yesterday on functions, `05-python-functions`, one exercise
asked you to write a function to transcribe DNA to RNA. An example of
that function is implemented in this directory in a file named
`transcribe.py`. In that lesson, there were two tests to check your
work:

```python
transcribe('ATGC') == 'UACG'
transcribe('ATGCAGTCAGTGCAGTCAGT') == 'UACGUCAGUCACGUCAGUCA'
```

Convert these to a proper test and place it the file `test_transcribe.py`.
Next, add a few tests of your own and run the test suite with nosetests.
