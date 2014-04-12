---
layout: lesson
root: ../..
title: Branch Out
---
####Learning Branching

#####Branches defined: 

* Branching is a way for you to create a new working copy of a project without affecting the main project.

#####Benefits of branching:

* You can change a project radically and save it, but always have the ability to switch back to your stable copy.
* Other people or collaborators can modify the project without affecting the main project. 

Let's explore how this tool can help with Dracula, Mummy and Wolfman's new home, `mars.txt`.

While `mars.txt` is great for Dracula and Mummy, Wolfman is having a hard time with all the changes. The two moons are really cramping his style! He decides to fix the problem by proposing `shelter.txt`

Follow along by typing the following into your terminal: 

<div class="in" markdown="1">
~~~
$ Touch shelter.txt

$ Nano shelter.txt
~~~
</div>

**Copy and Paste the following into shelter.txt:**

<div class="in" markdown="1">
~~~
Bedrooms:
3

* For Dracula: Coffin 10x10
* For Mummy: Sarcaphogus 10x10 
* For me: Bedroom with no windows 80x80
~~~
</div>

Let's check that our file is there.

<div class="in" markdown="1">
~~~
ls
~~~
</div>

####Great!
Wolfman isn't sure these dimensions will be agreeable, and Mummy and Dracula aren't the best people to upset when you're the only warm-blooded creature on Mars! So instead of pushing his changes to the main project, he can make his changes there to ask Dracula and Mummy for feedback before he without disrupting their other work.

<div class="in" markdown="1">
~~~
git branch shelter

git add shelter.txt
git status
git commit shelter.txt
git status
git push
~~~
</div>

We just pushed `shelter.txt` to a new branch of the main project! Let's see if that worked!

<div class="in" markdown="1">
~~~
git branch
~~~
</div>

We see that there are two branches possible:
* shelter
  master

The master is our main project! We're not there so any changes we push won't be seen by Dracula or Mummy! Phew! That's a relief. 

Let's make some changes that will be more fair to everyone:
<div class="in" markdown="1">
~~~
nano shelter.txt
~~~
</div>

Edit the numbers in shelter.txt to be more fair as so:
<div class="in" markdown="1">
~~~
Bedrooms:
3

* For Dracula: Coffin 33.3x33.3
* For Mummy: Sarcaphogus 33.3x33.3 
* For me: Bedroom with no windows 33.3x33.3
~~~
</div>

Ah, that's much more fair. Wolfman is now ready to commit these changes to the main project. So he needs to swtich branches to do so by using the command `git checkout`.

Let's try that on for size:

<div class="in" markdown="1">
~~~
git checkout master
git branch
~~~
</div>

We should see that we are back to the main project, known as the master in our output, as shown below:

* master
  shelter

Great! Now we can "push our changes to the master branch!"


<div class="in" markdown="1">
~~~
git add shelter.txt
git status
git commit sheter.txt
git status
git push
~~~
</div>

Let's checkout our repository online at github.com to ensure we have added the right dimensions to shelter.txt!

