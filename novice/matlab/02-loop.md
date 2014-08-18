---
layout: lesson
root: ../..
---

<div class="Objectives">

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

## Saving Our Work

So far, we've been entering commands one-by-one to get
MATLAB to do things for us. But what if we want to repeat 
our analysis?

Sure, it's only a handful of commands, and it shouldn't take
us more than a few minutes. But we'll quickly find ourselves
doing more complex analyses, and we'll need our results to
be more easily reproducible.

### Writing Scipts

In addition to running MATLAB commands one-by-one, we can 
also write several commands in a _script_. A MATLAB scipt
is just a text file with a `.m` extension. We've written 
commands to load data from a `.csv` file and
dsiplays some statisitics about that data. Let's 
put those commands in a script:

~~~
% script analyze.m

patient_data = csvread('inflammation-01.csv');

disp([Analyzing "inflammation-01.csv: "])
disp(['Maximum inflammation: ', num2str(max(patient_data(:)))]);
disp(['Minimum inflammation: ', num2str(min(patient_data(:)))]);
disp(['Standard deviation: ', num2str(std(patient_data(:)))]);
~~~

This script needs to be visible to MATLAB. MATLAB keeps an eye
on several locations on your [filesystem](../../gloss.html#filesystem).
The most convenient of these places to save
a script is generally in the _current directory_.
To find out the current directory, use the `pwd` command:

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

We've also written commands to create plots. MATLAB let's us save those as
images on disk.

~~~
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

% save plot to disk as png image:
print -dpng "inflammation-01-analysis.png"
~~~

{:class="in"}

Let's extend our `analyze` script with these commands.

~~~
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

% save plot to disk as png image:
print -dpng "inflammation-01-analysis.png"
~~~

## Analyzing Multiple Datasets

We have a dozen data sets right now, and more on the way. We want
to create plots for all our data sets without repeating the
above commands each time. To do that we'll have to learn how to
get the computer to repeat things.


### For loops

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
for i = 1:4
    word(i)
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
the name `i`.

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

for i in 3
    disp(word(i));
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

for i = 1:length(word)
    disp(word(i));
end
~~~
{:class="in"}

~~
a
l
u
m
i
n
i
u
m
~~
{:class="out"}

This is much more robust code, as it can deal indentically with
words of arbitrary length. Here's another loop that repeatedly 
updates the variable `len`.:

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
executed five times. The first time around, `len` is zero (the
value assigned to it before the loop. The loop body adds 1 to the old
value of `len`, producing 1, and updates `len` to refer to that new
value.
The next time around, `vowel` is `e`, and `len` is `, so `len` is
updated to 2.
After three more updates, `len` is 5; since thre's nothing left in
`aeiou` for MATLAB to process, the loop finishes and the
`disp` statement tells us our final answer.

Note that a loop variable is just a variable that's being used to record
progrss in aloop. It still exists after the loop is over, and we can re-use
variables previously defined as loop variables as well:

~~~
disp(vowel)
~~~
{:class="in"}

~~~
u
~~~
{:class="out"}

After the loop, `vowel` refers to the last value in `'aeiou'`, i.e., `'u'`.

## Processing Multiple Files

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

This is close, but not quite right. The `sprintf` function works a
lot like the `printf` function that you might have used if you've
programmed in C. The `sprintf` function is useful when we want to
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

it substitutes the `%d` in the template `'inflammation-%d.csv'`, with the
value of `i`, i.e., 1. The resulting string is `'inflammation-1.csv'`,
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

Key Points

* Write _scripts_ that store related MATLAB commands, and make
their results more reproducible.

* Save images generated by MATLAB using the `print` function.

* Use `for variable = collection` to process the elements
  of a collection (a MATLAB array), one at a time.

* Use the `length` command to determine the length of a MATLAB array.

* Use the `sprintf` function to generate a string based on a template.

