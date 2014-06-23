---
layout: lesson
root: ../..
title: Numbers
---
Let's start by looking at how numbers are stored.
If we only have the two digits 0 and 1,
the natural way to store a positive integer is to use base 2,
so 1001<sub>2</sub> is
(1&times;2<sup>3</sup>)+(0&times;2<sup>2</sup>)+(0&times;2<sup>1</sup>)+(1&times;2<sup>0</sup>) = 9<sub>10</sub>.
It's equally natural to extend this scheme to negative numbers by reserving one bit for the sign.
If, for example, we use 0 for positive numbers and 1 for those that are negative,
+9<sub>10</sub> would be 01001<sub>2</sub> and -9<sub>10</sub> would be 11001<sub>10</sub>.

There are two problems with this.
The first is that this scheme gives us two representations for zero (00000<sub>2</sub> and 10000<sub>2</sub>).
This isn't necessarily fatal,
but any claims this scheme has to being "natural" disappear when we have to write code like:

~~~
if (length != +0) and (length != -0)
~~~
{:class="in"}

As for the other problem,
it turns out that the circuits needed to do addition and other arithmetic on this
[sign and magnitude representation](../../gloss.html#sign-and-magnitude)
are more complicated than the hardware needed for another called
[two's complement](../../gloss.html#twos-complement).
Instead of mirroring positive values,
two's complement rolls over when going below zero,
just like a car's odometer.
If we're using four bits per number,
so that 0<sub>10</sub> is 0000<sub>2</sub>,
then -1<sub>10</sub> is 1111<sub>2</sub>.
-2<sub>10</sub> is 1110<sub>2</sub>,
-3<sub>10</sub> is 1101<sub>2</sub>,
and so on until we reach the most negative number we can represent,
1000<sub>2</sub>, which is -8.
Our representation then wraps around again, so that 0111<sub>2</sub> is 7<sub>10</sub>.

This scheme isn't intuitive,
but it solves sign and magnitude's "double zero" problem,
and the hardware to handle it is faster and cheaper.
As a bonus,
we can still tell whether a number is positive or negative by looking at the first bit:
negative numbers have a 1, positives have a 0.
The only odd thing is its asymmetry:
because 0 counts as a positive number,
numbers go from -8 to 7, or -16 to 15, and so on.
As a result, even if `x` is a valid number, `-x` may not be.

Finding a good representation for real numbers
(called [floating point numbers](../../gloss.html#float-point-number),
since the decimal point can move around)
is a much harder problem.
The root of the problem is that
we cannot represent an infinite number of real values with a finite set of bit patterns.
And unlike integers,
no matter what values we *do* represent,
there will be an infinite number of values between each of them that we can't.

Floating point numbers are usually represented using sign, magnitude, and an exponent.
In a 32-bit word,
the IEEE 754 standard calls for 1 bit of sign,
23 bits for the magnitude (or *mantissa*),
and 8 bits for the exponent.
To illustrate the problems with floating point,
we'll use a much dumber representation:
we'll only worry about positive values without fractional parts,
and we'll only use 3 for the magnitude and 2 for the exponent.

<!--- Remove this style when vertical headers was supported by pandoc:
https://github.com/jgm/pandoc/issues/1359 -->
<style>
.table-exponent td {
    width:17%;
}
.table-exponent td.table-exponent-header {
    font-weight: bold;
}
</style>

<!--- Merge cells around "Exponent" when colspan was supported by pandoc:
https://github.com/jgm/pandoc/issues/1340 -->
<table class="table table-striped table-exponent">
<tr><td></td>        <td>   </td><td></td><td class="table-exponent-header">Exponent</td><td></td><td></td></tr>
<tr><td></td>        <td>   </td><td class="table-exponent-header">00</td><td class="table-exponent-header">01</td><td class="table-exponent-header">10</td><td class="table-exponent-header">11</td></tr>
<tr><td></td>        <td class="table-exponent-header">000</td><td> 0</td><td> 0</td><td> 0</td><td> 0</td></tr>
<tr><td></td>        <td class="table-exponent-header">001</td><td> 1</td><td> 2</td><td> 4</td><td> 8</td></tr>
<tr><td></td>        <td class="table-exponent-header">010</td><td> 2</td><td> 4</td><td> 8</td><td>16</td></tr>
<tr><td class="table-exponent-header">Mantissa</td><td class="table-exponent-header">011</td><td> 3</td><td> 6</td><td>12</td><td>24</td></tr>
<tr><td></td>        <td class="table-exponent-header">100</td><td> 4</td><td> 8</td><td>16</td><td>32</td></tr>
<tr><td></td>        <td class="table-exponent-header">101</td><td> 5</td><td>10</td><td>20</td><td>40</td></tr>
<tr><td></td>        <td class="table-exponent-header">110</td><td> 6</td><td>12</td><td>24</td><td>48</td></tr>
<tr><td></td>        <td class="table-exponent-header">111</td><td> 7</td><td>14</td><td>28</td><td>56</td></tr>
</table>

The table above
shows the values that we can represent this way.
Each one is the mantissa times two to the exponent.
For example, the decimal values 48 is binary 110 times 2 to the binary 11 power,
which is 6 times 2 to the third,
or 6 times 8.
(Note that real floating point representations like the IEEE 754 standard
don't have the redundancy shown in this table,
but that doesn't affect our argument.)

The first thing you should notice is that there are a lot of values we *can't* store.
We can do 8 and 10, for example, but not 9.
This is exactly like the problems hand calculators have with fractions like 1/3:
in decimal, we have to round that to 0.3333 or 0.3334.

But if this scheme has no representation for 9,
then 8+1 must be stored as either 8 or 10.
This raises an interesting question:
if 8+1 is 8, what is 8+1+1?
If we add from the left, 8+1 is 8, plus another 1 is 8 again.
If we add from the right, though, 1+1 is 2, and 2+8 is 10.
Changing the order of operations can make the difference between right and wrong.
There's no randomness involved&mdash;a particular order of operations
will always produce the same result&mdash;but
as the number of steps increases,
so too does the difficulty of figuring out what the best order is.

This is the sort of problem that numerical analysts spend their time on.
In this case, if we sort the values we're adding, then add from smallest to largest,
it gives us a better chance of getting the best possible answer.
In other situations,
like inverting a matrix,
the rules are much more complicated.

Here's another observation about our uneven number line:
the spacing between the values we can represent is uneven,
but the relative spacing between each set of values stays the same,
i.e., the first group is separated by 1, then the separation becomes 2, then 4, then 8,
so that the ratio of the spacing to the values stays roughly constant.
This happens because we're multiplying the same fixed set of mantissas by ever-larger exponents,
and it points us at a couple of useful definitions.

The [absolute error](../../gloss.html#absolute-error) in some approximation
is simply the absolute value of the difference between the actual value and the approximation.
The [relative error](../../gloss.html#relative-error),
on the other hand,
is the ratio of the absolute error to the value we're approximating.
For example, if we're off by 1 in approximating 8+1 and 56+1,
the absolute error is the same in both cases,
but the relative error in the first case is 1/9 = 11%,
while the relative error in the second case is only 1/57 = 1.7%.
When we're thinking about floating point numbers,
relative error is almost always more useful than absolute error.
After all,
it makes little sense to say that we're off by a hundredth when the value in question is a billionth.

To see why this matters, let's have a look at a little program:

~~~
nines = []
sums = []
current = 0.0
for i in range(1, 10):
    num = 9.0 / (10.0 ** i)
    nines.append(num)
    current += num
    sums.append(current)

for i in range(len(nines)):
    print '%.18f %.18f' % (nines[i], sums[i])
~~~
{:class="in"}

The loop runs over the integers from 1 to 9 inclusive.
Using those values, we create the numbers 0.9, 0.09, 0.009, and so on, and put them in the list `vals`.
We then calculate the sum of those numbers.
Clearly, this should be 0.9, 0.99, 0.999, and so on.
But is it?

<table class="table table-striped">
<tr><td>1</td><td>0.900000000000000022</td><td>0.900000000000000022</td></tr>
<tr><td>2</td><td>0.089999999999999997</td><td>0.989999999999999991</td></tr>
<tr><td>3</td><td>0.008999999999999999</td><td>0.998999999999999999</td></tr>
<tr><td>4</td><td>0.000900000000000000</td><td>0.999900000000000011</td></tr>
<tr><td>5</td><td>0.000090000000000000</td><td>0.999990000000000046</td></tr>
<tr><td>6</td><td>0.000009000000000000</td><td>0.999999000000000082</td></tr>
<tr><td>7</td><td>0.000000900000000000</td><td>0.999999900000000053</td></tr>
<tr><td>8</td><td>0.000000090000000000</td><td>0.999999990000000061</td></tr>
<tr><td>9</td><td>0.000000009000000000</td><td>0.999999999000000028</td></tr>
</table>

Here are our answers.
The first column is the loop index;
the second, what we actually got when we tried to calculate 0.9, 0.09, and so on,
and the third is the cumulative sum.

The first thing you should notice is that the very first value contributing to our sum is already slightly off.
Even with 23 bits for a mantissa,
we cannot exactly represent 0.9 in base 2,
any more than we can exactly represent 1/3 in base 10.
Doubling the size of the mantissa would reduce the error,
but we can't ever eliminate it.

The second thing to notice is that our approximation to 0.0009 actually appears accurate,
as do all of the approximations after that.
This may be misleading, though:
after all,
we've only printed things out to 18 decimal places.
As for the errors in the last few digits of the sums,
there doesn't appear to be any regular pattern in the way they increase and decrease.

This phenomenon is one of the things that makes testing scientific programs hard.
If a function uses floating point numbers,
what do we compare its result to
if we want to check that it's working correctly?
If we compared the sum of the first few numbers in `vals` to what it's supposed to be,
the answer could be `False`,
even if we're initializing the list with the right values,
and calculating the sum correctly.
This is a genuinely hard problem,
and no one has a good generic answer.
The root of our problem is that we're using approximations,
and each approximation has to be judged on its own merits.

There are things you can do, though.
The first rule is,
compare what you get to analytic solutions whenever you can.
For example,
if you're looking at the behavior of drops of liquid helium,
start by checking your program's output on a stationary spherical drop in zero gravity.
You should be able to calculate the right answer in that case,
and if your program doesn't work for that,
it probably won't work for anything else.

The second rule is to compare more complex versions of your code to simpler ones.
If you're about to replace a simple algorithm for calculating heat transfer with one that's more complex,
but hopefully faster,
don't throw the old code away.
Instead,
use its output as a check on the correctness of the new code.
And if you bump into someone at a conference who has a program that can calculate some of the same results as yours,
swap data sets:
it'll help you both.

The third rule is, never use `==` (or `!=`) on floating point numbers,
because two numbers calculated in different ways will probably not have exactly the same bits.
Instead,
check to see whether two values are within some tolerance,
and if they are,
treat them as equal.
Doing this forces you to make your tolerances explicit,
which is useful in its own right
(just as putting error bars on experimental results is useful).
