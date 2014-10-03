---
layout: lesson
root: ../..
title: Reusing Code
---


We've written MATLAB commands to compute statistics about our
data and generate some plots to visualize the results. We're now
faced with the following problems:

**Problem 1**

So far, we've typed in commands one-by-one on the command line
to get MATLAB to do things for us. But what if we want to repeat
our analysis? Sure, it's only a handful of commands,
and typing them in shouldn't take
us more than a few minutes. But if we forget a step or make a mistake,
we'll waste time rewriting commands. Also, we'll quickly find ourselves
doing more complex analyses, and we'll need our results to
be more easily reproducible.

**Problem 2**

We also have to do this analysis for every one of our dozen
datasets. And we need a better way than typing out commands
for each one, because we'll find ourselves writing a lot
of duplicate code. Remember, code that is repeated in
two or more places will eventually be wrong in at least one.
If we  make changes in the way we analyze our datasets,
we have to introduce that change in every copy of our code.

There's a common theme in the two problems presented above---duplicate code.
In problem 1, we're rewriting code every time we want
to perform the same analysis several times. In problem 2,
we're rewriting code
every time we want to perform several similar analyses.
To avoid writing all this duplicate code, we have to teach MATLAB to

* Remember our commands

* Repeat our commands

<div class="objectives" markdown="1">

### Objectives

    * Learn how to write and save MATLAB scripts.
    * Learn how to save MATLAB plots to disk.
    * Explain what a for loop does.
    * Correctly write for loops that repeat simple commands.
    * Trace changes to a loop variable as the loop runs.
    * Use a for loop to process multiple files.

</div>

<!-- FIXME: The 0 lesson should talk about the colon operator -->
<!-- FIXME: Need to mention that strings are arrays -->
<!-- FIXME: Decide on consistent string---single or double quotes? -->

### Saving Our Work

### Writing Scripts

In addition to running MATLAB commands one-by-one on the
command line, we can
also write several commands in a _script_. A MATLAB script
is just a text file with a `.m` extension. We've written
commands to load data from a `.csv` file and
displays some statistics about that data. Let's
put those commands in a script called `analyze.m`:

~~~
% script analyze.m

patient_data = csvread('inflammation-01.csv');

disp([Analyzing "inflammation-01.csv: "])
disp(['Maximum inflammation: ', num2str(max(patient_data(:)))]);
disp(['Minimum inflammation: ', num2str(min(patient_data(:)))]);
disp(['Standard deviation: ', num2str(std(patient_data(:)))]);
~~~

Before we can use it, we need to make sure that this file is
visible to MATLAB. MATLAB doesn't know about all the files on your
computer, but it keeps an eye on several directories.
The most convenient of these directories is generally the
"working directory", or "current directory". To find out the
working directory, use the `pwd` command:

~~~
pwd
~~~
{:class="in"}

As you might have guessed, `pwd` stands for "print working directory".

Once you have a script saved in a location that MATLAB knows about,
you can get MATLAB to run those commands by typing in the name
of the script (without the `.m`) in the MATLAB command line:

~~~
analyze
~~~
{:class="in"}

~~~
Maximum inflammation: 20
Minimum inflammation: 0
Standard deviation: 4.7219
~~~
{:class="out"}

### Saving Images

We've also written commands to create plots:

~~~
ave_inflammation = mean(patient_data, 1);

plot(ave_inflammation);
ylabel("average")
~~~
{:class="in"}


MATLAB let's us save those as
images on disk:

~~~
% save plot to disk as png image:
print -dpng "average.png"
~~~
{:class="in"}

Let's extend our `analyze` script with commands to
create and save plots:

~~~
% script analyze.m

patient_data = csvread('inflammation-01.csv');

disp(['Maximum inflammation: ', num2str(max(patient_data(:)))]);
disp(['Minimum inflammation: ', num2str(min(patient_data(:)))]);
disp(['Standard deviation: ', num2str(std(patient_data(:)))]);

ave_inflammation = mean(patient_data, 1);

subplot(1, 3, 1);
plot(ave_inflammation);
ylabel("average")

subplot(1, 3, 2);
plot(max(patient_data, [], 1));
ylabel("max")

subplot(1, 3, 3);
plot(min(patient_data, [], 1));
ylabel("min")

% save plot to disk as svg image:
print -dpng "patient_data-01.png"
~~~


#### The Colon Operator
<!-- FIXME: Maybe this should be in a box? -->

You can use the `:` (colon) operator to generate
sequences in MATLAB:

~~~
4:10
~~~
{:class="in"}

~~~
ans =

    4    5    6    7    8    9   10
~~~
{:class="out"}


~~~
2.5:0.25:5
~~~
{:class="in"}

~~~
ans =

 Columns 1 through 8:

    2.5000    2.7500    3.0000    3.2500    3.5000    3.7500    4.0000    4.2500

 Columns 9 through 11:

    4.5000    4.7500    5.0000

~~~
{:class="out"}

### Analyzing Multiple Datasets

We have a dozen data sets right now, and more on the way. We want
to create plots for all our data sets without repeating the
above commands each time. To do that we'll have to learn how to
get the computer to repeat things.


### for loops

Suppose we want to print each character in the word "lead" on
a line of its own. One way is to use four `disp` statements:

~~~
word = 'lead';

disp(word(1));
disp(word(2));
disp(word(3));
disp(word(4));
~~~
{:class="in"}

~~~
l
e
a
d
~~~
{:class="out"}

But that's a bad approach for two reasons:

1. It doesn't scale: if we want to print the characters
in a string that's hundreds of letters long, we'd be better
off typing them in.

2. It's fragile: if we change `word` to a longer string,
it only
prints part of the data, and if we change it to
a shorter one,
it produces an error, because we're asking for characters
that don't exist.

~~~
word = 'tin'

disp(word(1));
disp(word(2));
disp(word(3));
disp(word(4));
~~~
{:class="in"}

~~~
error: A(I): index out of bounds; value 4 out of bound 3
~~~
{:class="out"}


There's a better approach:

~~~
for letter = 1:4
    disp(word(letter))
end
~~~
{:class="in"}

~~~
l
e
a
d
~~~
{:class="out"}

This improved version uses a [for loop](../../gloss.html#for-loop) to
repeat an operation---in this case, printing---once for each element in
an array.

The general form of a for loop is:

~~~
for variable = collection
    do things with variable
end
~~~

The for loop executes the commands in the
[loop body](../../gloss.html#loop-body)
for every value in the array `collection`. This
value is called the [loop variable](../../gloss.html#loop-variable),
and we can call it whatever we like. In our example, we gave it
the name `letter`.

We have to terminate the loop body with the `end` keyword, and
we can have as many commands as we like in the loop body. But we
have to remember that they will all
be repeated as many times as
there are values in `collection`.

Our for loop has made our code more scalable, and less fragile.
There's still one little thing about it that should bother us.
For our loop to deal appropriately with shorter or longer words,
we have to change the first line of our loop by hand:

~~~
word = 'tin';

for letter = 1:3
    disp(word(letter));
end
~~~
{:class="in"}

~~~
t
i
n
~~~
{:class="out"}


Although this works, it's not the best way to write our loop:

* We might update `word` and forget to modify the loop to reflect that
  change.

* We might make a mistake while counting the number of letters in
  `word`.

Fortunately, MATLAB provides us with a convenient function to
write a better loop:

~~~
word = 'aluminium'

for letter = 1:length(word)
    disp(word(letter));
end
~~~
{:class="in"}

~~~
a
l
u
m
i
n
i
u
m
~~~
{:class="out"}

This is much more robust code, as it can deal indentically with
words of arbitrary length. Here's another loop that repeatedly
updates the variable `len`:

~~~
len = 0
for vowel = 'aeiou'
    len = len + 1;
end

disp(['Number of vowels: ', num2str(len)])
~~~
{:class="in"}


It's worth tracing the execution of this little program step by step.
Since there are five characters in "aeiou", the loop body will be
executed five times. When we enter the loop, `len` is zero - the
value assigned to it beforehand. The first time through, the loop body adds 1 to the old
value of `len`, producing 1, and updates `len` to refer to that new
value.
The next time around, `vowel` is `e`, and `len` is 1, so `len` is
updated to 2.
After three more updates, `len` is 5; since there's nothing left in
`aeiou` for MATLAB to process, the loop finishes and the
`disp` statement tells us our final answer.

Note that a loop variable is just a variable that's being used to record
progress in a loop. It still exists after the loop is over, and we can re-use
variables previously defined as loop variables as well:

~~~
disp(vowel)
~~~
{:class="in"}

~~~
u
~~~
{:class="out"}

After the loop, `vowel` refers to the last value in `aeiou`, i.e., `u`.

#### Challenges

1. Write a loop that spells the word "aluminum," adding one letter at a time: 
 
   ~~~
   a
   al
   alu
   alum
   alumi
   alumin
   aluminu
   aluminum
   ~~~
   {:class="out"}

1. Matlab uses the caret (`^`) to perform exponentiation:

   ~~~
   disp(5^3)
   ~~~
   {:class="in"}
   ~~~
   125
   ~~~
   {:class="out"}

   Let `b` be the base of the number and `x` the exponent. Write a loop to compute `b^x`. Check your result for `b = 4` and `x = 5`.

1. In Matlab, the colon operator (`:`) accepts a [stride](../../gloss.html#stride) or skip argument between the start and stop:

   ~~~
   disp(1:3:11)
   ~~~
   {:class="in"}
   ~~~
   1 4 7 10
   ~~~
   {:class="in"}
   ~~~
   disp(11:-3:1)
   ~~~
   {:class="in"}
   ~~~
   11 8 5 2
   ~~~
   {:class="out"}

   Using this, write a loop to print the letters of "aluminum" in reverse order, one letter per line.

   ~~~
   m
   u
   n
   i
   m
   u
   l
   a
   ~~~ 
   {:class="out"}

**Extra Challenge**: Reverse the string `abracadabra` without a loop, using only indexing and the colon operator.

### Processing Multiple Files

We now have almost everything we need to process
multiple data files with our `analyze` script. You'll notice that our
data files are named `inflammation-01.csv`, `inflammation-02.csv`, etc.
Let's write a loop that tries to print the names of each one of
our files:

~~~
for i = 1:12
    file_name = sprintf('inflammation-%d.csv', i);
    disp(file_name);
end
~~~
{:class="in"}

~~~
inflammation-1.csv
inflammation-2.csv
inflammation-3.csv
inflammation-4.csv
inflammation-5.csv
inflammation-6.csv
inflammation-7.csv
inflammation-8.csv
inflammation-9.csv
inflammation-10.csv
inflammation-11.csv
inflammation-12.csv
~~~
{:class="out"}

This is close, but not quite right.
The `sprintf` function is useful when we want to
generate MATLAB strings based on a _template_. In our case,
that template is the string `inflammation-%d.csv`. `sprintf`
generates a new string from our template by replacing the `%d` with
the data referred to by our second argument, `i`.

Again, let's trace the execution of our loop: in the beginning of our
loop, `i` starts by referring to
the value 1. So, when MATLAB executes the command

~~~
file_name = sprintf('inflammation-%d.csv', i);
~~~

it substitutes the `%d` in the template `inflammation-%d.csv`, with the
value of `i`, i.e., 1. The resulting string is `inflammation-1.csv`,
which is assigned to the variable `file_name`. The `disp` command
prints that string to screen. The second time around, `sprintf`
generates the string `inflammation-2.csv`, which is assigned to the
variable `file_name`, and printed to screen. And
so on, till `i` finally refers to the value 12.

Remember that there's a mistake. Our files are actually
named `inflammation-01.csv`, `inflammation-02.csv`, etc. To get it
right, we have to modify our template:

~~~
for i = 1:12
    file_name = sprintf('inflammation-%02d.csv', i);
    disp(file_name);
end
~~~
{:class="in"}

~~~
inflammation-01.csv
inflammation-02.csv
inflammation-03.csv
inflammation-04.csv
inflammation-05.csv
inflammation-06.csv
inflammation-07.csv
inflammation-08.csv
inflammation-09.csv
inflammation-10.csv
inflammation-11.csv
inflammation-12.csv
~~~
{:class="out"}

We've replaced `%d` in our earlier template with `%02d`. With this,
we're specifying that we want our data to be displayed with a minimum
width of 2 characters, and that we want to _pad_ with 0 for data that
isn't at least 2 digits long.

We're now ready to modify `analyze.m` to process multiple data files:

~~~
% script analyze.m

for i = 1:3

    % Generate strings for file and image names:
    file_name = sprintf('inflammation-%02d.csv', i);
    img_name = sprintf ('patient_data-%02d.svg', i);

    patient_data = csvread(file_name);
    ave_inflammation = mean(patient_data, 1);

    figure()

    subplot(1, 3, 1);
    plot(ave_inflammation);
    ylabel('average')

    subplot(1, 3, 2);
    plot(max(patient_data, [], 1));
    ylabel('max')

    subplot(1, 3, 3);
    plot(min(patient_data, [], 1));
    ylabel('min')

    print(img_name);
    close();

end
~~~

Remember that to run our script, we simply type in its name in the
command line:

~~~
analyze
~~~
{:class="in"}

<!-- FIXME: put the plots here -->

Sure enough, the maxima of these data sets show exactly the same ramp
as the first, and their minima show the same staircase structure.

#### Challenges

1. In cases where our file names don't follow such a regular pattern, we might
   want to process all files that end with a given extension, say `.csv`. At the
   command line we could get this list of files by using a
   [wildcard](../../gloss.html#wildcard):

   ~~~
   ls *.csv
   ~~~
   {:class="in"}

   Thankfully, Matlab also has `ls`, though it returns a single long string:
   
   ~~~
   filestr = ls('*.csv')
   ~~~
   {:class="in"}
   
   ~~~
   inflammation-01.csv inflammation-04.csv inflammation-07.csv inflammation-10.csv
   inflammation-02.csv inflammation-05.csv inflammation-08.csv inflammation-11.csv
   inflammation-03.csv inflammation-06.csv inflammation-09.csv inflammation-12.csv
   ~~~
   {:class="out"}
   
   To turn this string into an array we can loop over (actually, a 
   [Cell Array](http://www.mathworks.com/help/matlab/cell-arrays.html)),
   we need to "split" the string at each occurrence of whitespace:
   
   ~~~
   file_string = ls('*.csv');
   file_list = strsplit(file_string)
   ~~~
   {:class="in"}
   
   ~~~
   file_list = 
   
     Columns 1 through 3
   
       'inflammation-01.csv'    'inflammation-04.csv'    'inflammation-07.csv'
   
     Columns 4 through 6
   
       'inflammation-10.csv'    'inflammation-02.csv'    'inflammation-05.csv'
   
     Columns 7 through 9
   
       'inflammation-08.csv'    'inflammation-11.csv'    'inflammation-03.csv'
   
     Columns 10 through 13
   
       'inflammation-06.csv'    'inflammation-09.csv'    'inflammation-12.csv'    ''
   ~~~
   {:class="out"}

   Using this trick, rewrite the `analyze` script to analyze all `csv` files in
   the current directory. Be careful of the empty string `''` at the end of
   `file_list`!

<div class="keypoints" markdown="1">
#### Key Points

* Write MATLAB scripts to reuse code and make your results reproducible.

* Save images generated by MATLAB using the `print` function.

* Use a for loop: `for variable = collection`, to process the elements
  of a collection (a MATLAB array), one at a time.

* Use the `:` (colon) operator to generate sequences.

* Use the `length` command to determine the length of a MATLAB array.

* Use the `sprintf` function to generate a string based on a template.
</div
