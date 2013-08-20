The following exercises do not contain solutions. Yet. Instead, we will be
asking you to submit your solutions to these exercises and then we will post up
solutions at the start of next week. We encourage you to discuss your approaches
or solutions on the course forum!

To submit your exercises, please create a `testing` folder in your personal
folder in the course repository. Place all of the code and files for theses
exercises in that folder and be sure to check it in.


## Exercise 1: Mileage

The function 'convert_mileage' converts miles per gallon (US style) toÂ liters
per 100 km (metric style):

```python
    gal_to_litre = 3.78541178
    mile_to_km = 1.609344
    
    def convert_mileage(mpg):
        '''Converts miles per gallon to liters per 100 km'''
        litres_per_100_km = 100 / mpg / mile_to_km * gal_to_litre
        return litres_per_100_km
```

Create a subdirectory in your version control directory called `testing`, then
copy this function into a file in that directoryÂ called `mileage.py`. Add more
code to that file to repeatedly ask the user for a mileage in miles per gallon,
and output the mileage in liters per 100 km, until the user enters the string
"`q`". You will need to use the `float()` function to convert from string to a
floating point number. Use the '`if __name__ == "__main__":`' trick to ensure
that the module can be imported without executing your testing code.

1. Copy `mileage.py` to create `tryexcept.py` Add a try/except block to the new
program to display a helpful message instead of crashing when users enter
invalid input (such as the number "0" or the name of their favorite hockey
team). 

2. Reading the function again, you realize that accepting 0 or negative values
make no sense and should be reported as an error. Look at the exceptions defined
in the `exceptions` module (use the built-in `help(...)` or `dir(...)`
functions) and decide which of Python's built-in exceptions is most appropriate
to use for invalid input. Create a copy of 'tryexcept.py' called 'raiser.py'
that raises this exception; modify the main body of your program to catch it;
and add a comment inside the file explaining why you chose the exception you
did. (Note: you have to call this file `raiser.py`, not `raise.py` because
'import raise' is an error. Â Can you see why?)

3. [According to
Google](http://www.google.ca/search?q=20+miles+per+gallon+in+litres+per+100+km&gbv=1),
20 miles per gallon are equivalent to 11.7607292 liters per 100 km. Use these
values to write a unit test. Keep in mind that these floating values are subject
to truncation and rounding errors. Save the test case in a file called
`test_mileage.py` and run it using the `nosetests` command. Â Note:
`test_mileage.py` should use '`from raiser import convert_mileage`' to get the
final version of your mileage conversion function.

4. Now add a second test case, for 40 miles per gallon equivalent to 5.88036458
liters per 100 km and run the tests again. Â Unless you have already fixed the
error that was present in the initial function, your test should fail. Â Find
and fix the error; submit your new function in a file called 'final_mileage.py'. 


## Exercise 2: Testing Averages

The results of a set of experiments are stored in a file, where the _i-th_ line
stores the results of the _i-th_ experiment as a comma-separated list of
integers. A student is assigned the task of finding the experiment with the
smallest average value. She writes the following code: 

```python
    def avg_line(line):
           values = line.split(',')
           count = 0
           total = 0
           for value in values:
                total += int(value)
                count += 1
           return total / count
    
       def min_avg(file_name):
           contents = open(file_name)
           averages = []
           for (i, line) in enumerate(contents):
               averages.append((avg_line(line), i))
           contents.close()
           averages.sort()
           min_avg, experiment_number = averages[0]
           return experiment_number
```

1. Refactor `min_avg` so that it can be tested without depending on external
files. Submit your code in a file called `first_averages.py`.

2. Write Nose test cases for both functions. Consider what should happen if the
file is empty. Submit your tests in a file called `test_first_averages.py`.
Note: you may assume for now that all input is well formatted, i.e., you do
_not_ have to worry about empty lines, lines containing the names of hockey
teams, etc.

3. The given specification is ambiguous: what should the result be if two or
more experiments are tied for the minimum average?  Copy 'first_averages.py' to
create a new file 'second_averages.py'; modify it to handle this case; add a
comment to the top explaining what rule you decided to use; and create a file
'test_second_averages.py' that tests your changes.

4. Another student proposed an alternative implementation of the min_avg
function:

```python
    def min_avg(file_name):
        contents = open(file_name).readlines()
        min_avg = avg_line(contents[0])
        min_index = 0
        for (i,line) in enumerate(contents):
            current_avg = avg_line(line)
            if current_avg <= min_avg:
                min_avg = current_avg
                min_index = i
        return min_index
```    

This implementation also finds an experiment with the smallest average, but
possibly a different one than the your function.  Modify your test cases so that
both your implementation and this one will pass.  (Hint: use the 'in' operator.)

5. One way to avoid the ambiguity of this specification is to define a
'min_avg_all' function instead, which returns a list with all the experiments
with the smallest average, and let the caller select one. Write tests for the
'min_avg_all' function, considering the following situations: an empty file,
exactly one experiment with minimum average, and more than one experiment with
minimum average. Keep in mind that in the last case, implementations could
return the list in different order. Write the tests the file "test_averages.py".
Use the same data as for the previous tests, if possible. You should use
variables to avoid code duplication. You don't need to implement the
'min_avg_all' function, but your test cases should be comprehensive enough to
serve as a specification for it.
