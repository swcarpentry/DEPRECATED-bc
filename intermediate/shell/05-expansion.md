---
layout: lesson
root: ../..
title: Filename Expansion
level: intermediate
---
<div class="objectives" markdown="1">
#### Objectives
*   FIXME
</div>

Now that we have our files renamed with the prefix "original-",
let's reverse course and take off the prefix.
By placing the variable's name within curly braces, e.g., `${filename}`,
you gain new powers&mdash;like the ability to modify the variable's value when you extract it.

<div class="in" markdown="1">
~~~
for filename in *.dat
do
    echo mv $filename ${filename#original-}
done
~~~
</div>

The `#` notation removes text from the beginning of a variable's value.
So this loop would print:

<div class="out" markdown="1">
~~~
mv original-basilisk.dat basilisk.dat
mv original-unicorn.dat unicorn.dat
~~~
</div>

And using `%` instead of `#` removes text from the end:

<div class="in" markdown="1">
~~~
for filename in *.dat
do
    echo mv $filename ${filename%.dat}
done
~~~
</div>

prints:

<div class="out" markdown="1">
~~~
mv original-basilisk.dat original-basilisk
mv original-unicorn.dat original-unicorn
~~~
</div>

> #### Avoid confusing variable names and text
> 
> Sometimes you may want to add something to the end of a variable's value.
> For example, you might add "backup" to the end of your files' names:
> 
> ~~~
> for filename in *.dat
> do
>     echo mv $filenamebackup
> done
> ~~~
> 
> Oops: each time through the loop the shell looks for a variable named `filenamebackup`, which doesn't exist.
> To avoid confusing our `$filename` variable with the text "backup",
> we can use curly braces:
>
> ~~~
> for filename in *.dat
> do
>     echo mv $filename ${filename}backup
> done
> ~~~
>
> which prints:
>
> ~~~
> mv original-basilisk.dat original-basilisk.datbackup
> mv original-unicorn.dat original-unicorn.datbackup
> ~~~

<div class="keypoints" markdown="1">
#### Key Points
*   FIXME
</div>
