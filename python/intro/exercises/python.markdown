The following exercises do not contain solutions.  Yet.  Instead, we will be
asking you to submit your solutions to these exercises and then we will post up
solutions at the start of next week.  We encourage you to discuss your
approaches or solutions on the course forum!

To submit your exercises, please create a `python` folder in your personal
folder in the course repository.  Place all of the code and files for theses
exercises in that folder and be sure to check it in. 

## Exercise 1: Short Problems


1. Write a bullet list description of what happens when Python evaluates the
statement `x += x - x` when `x` starts with the value 3.  Submit your answer in
a file called `evaluation.txt`.

2. In the United States, a car's fuel efficiency is measured in miles per
gallon.  In the metric system, it is usually measured in liters per 100
kilometers.
	
  1. Write a function called `convert_mileage` that converts from miles per
  gallon to liters per 100 kilometers.

	
  2. Test that your functions returns the right values for 20 and 40 miles per
  gallon.

	
  3. How did you figure out what the right value was?  How closely do the
  computer's results match the ones you expected?


  Submit your answer in a file called `mileage.py` (for part 1) and
  `mileage.txt` (for part 3).

3. The expression `1+'a'` is not legal in Python because 1 is a number and 'a'
is a string.  Do you think `1<'a'` is legal or not?  If it is legal, what is its
value, and why?  Submit your answer in a file called `comparison.txt`.

4. What is the value of each of the following expressions, and why? Submit your
answer in a file called `boolean.txt`.
  1. `True and not False`
  2. `True or True and False`
  3. `not True or not False`
  4. `True and not 0`
  5. `52 < 52.3`
  6. `1 + 52 < 52.3`
  7. `4 != 4.0`

5. Here is a possible definition of how `and` works:

   > If both operands are true, `and`'s result is its second value. If either is
   > false, `and`'s result is `False`.

   Is this the rule that Python actually uses?  If not, provide a counter example.
   Submit your answer in a file called `and.txt`.

6. You want an automatic wildlife camera to switch on if the light level is less
than 0.01 or if the temperature is above freezing, but not if both conditions
are true.  Your first attempt to write this is as follows:

   ```python    
       if (light < 0.01) or (temperature > 0.0):
          if (light < 0.01) and (temperature > 0.0):
              pass
          else:
              camera.on()
   ```

   A friend says that you could write this more simply as follows:

   ```python    
       if (light < 0.01) != (temperature > 0.0):
          camera.on()
   ```

   Is your friend right?  If so, explain why.  If not, give values for `light` and
   `temperature` that will produce different results for the two fragments of code.
   Either way, which do you find easier to understand, and why?  Submit your answer
   in a file called `wildlife.txt`.


## Exercise 2: Viewing Grades

A post-doc is teaching a biology course for the first time.  The university's
learning management system is "temporarily unavailable" while the system
administrators try to undo a recent upgrade, so the post-doc has to manage marks
herself for the foreseeable future.  She needs a program to help her see:
  * a list of all students with their marks;
  * the highest mark;
  * the lowest mark; and
  * the average mark.

Grades for each assignment are stored in a separate file.  Each file is
formatted like this:
    aturing: 86
    inewton: 67.5
    cdarwin: 90
    fnightingale: 99
    cvraman: 10

Write four functions that work as follows:
  1. `grade_show('ex3.dat', 'inewton')`: print inewton's grade in `ex3.dat`
  2. `grade_highest('ex3.dat')`: print the highest grade in `ex3.dat`
  3. `grade_lowest('ex3.dat')`: print the lowest grade in `ex3.dat`
  4. `grade_average('ex3.dat')`: print the average grade in `ex3.dat`

You may assume that all students' IDs are unique.  Submit your solution in a
file called `grade1.py`.  Please also submit a plain text file called
`grade1.txt` that explains ambiguities or omissions there are in this
specification, i.e., what do you have to implement that someone else might
reasonably implement differently?

## Exercise 3: Managing Grades

The post-doc also needs to be able to:
  * add a new student with a grade;
  * change a student's mark; and
  * remove a student

Create and submit another file called `grade2.py` containing functions that work
as follows:
  1. `grade_create('ex3.dat')`: create an empty grades file
  2. `grade_add('ex3.dat', 'inewton', 58)`: set inewton's grade in `ex3.dat` to
     58
  3. `grade_change('ex3.dat', 'inewton', 88)`: change inewton's grade in
     `ex3.dat` to 88
  4. `grade_remove('ex3.dat', 'inewton')`: remove inewton's grade from `ex3.dat`


## Exercise 4: Steganography

_Steganography_ is the art or science of hiding secret messages in plain sight.
For example, the first letter of every fifth line of the original edition of The
Tempest might spell out "Francis Bacon writ Shakesprs plays", or we could
introduce deliberate misspellings into an email message such that every wrong
letter, taken together, spelled out "Higgs Boson found".

Write a program that takes a short message and a positive integer N as
arguments, reads in a text file, and writes out the text with every N'th letter
changed to encode the message.  For example, suppose the file `example.txt`
contains:

```
Since 1997, Software Carpentry has taught scientists and engineers the concepts,
skills, and tools they need to use and build software more productively. All of
the content is freely available under a Creative Commons license, and we are
constantly adding and updating lectures, videos, and exercises.
```

If the program is run as:
> `python steg.py "Hacker Within" 5 example.txt`

then the output should be:

```
Hincea1997c SofkwareeCarprntry has Waughiscitntishs ani engnneers the concepts,
skills, and tools they need to use and build software more productively. All of
the content is freely available under a Creative Commons license, and we are
constantly adding and updating lectures, videos, and exercises.
```

Submit your solution in a file called `steg.py`.  Please also create a plain
text file called `steg.txt` containing answers to the following questions:
	
  1. What ambiguities or omissions are there in this specification? I.e., what
  do you have to implement that someone else might reasonably implement
  differently?

  2. This program's output is clearly garbled.  One way to make that less
  obvious is to increase the spacing between changes.  What are some others?

	
  3. How did you test your program?  How confident are you that it works
  correctly?  How confident do you think someone else could be that it works
  correctly?
