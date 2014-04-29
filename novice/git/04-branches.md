---
layout: lesson
root: ../..
title: Learning to Branch
---

#####Objectives:
*	Explain what branches and why they are helpful for individual and collaborative workflows.
*	Explain how to make a branch and navigate between branches.
*	Explain how to isolate work in a branch and merge it back to your main branch.


#####Branches defined: 

* Branching is a way for you to create a new working copy of a project without affecting the main project.

#####Benefits of branching:

* You can change a project radically and save it, but always have the ability to switch back to your stable copy. This can help you work with the project as you need solo.
* Other people or collaborators can modify the project without affecting the main project. This serves as an advantage for the group.

Let's explore how this tool can help with Dracula, Mummy and Wolfman's new home, `mars.txt`.

While `mars.txt` is great for Dracula and Mummy, Wolfman is having a hard time with all the changes. The two moons are really cramping his style! He decides to fix the problem by proposing `shelter.txt`

Follow along by typing the following into your terminal: 

<div class="in" markdown="1">
~~~
$ touch shelter.txt
$ nano shelter.txt
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

Wolfman isn't sure these dimensions will be agreeable, and Mummy and Dracula aren't the best people to upset when you're the only warm-blooded creature on Mars! So instead of pushing his changes to the main project, he can make his changes on a new "branch" to ask Dracula and Mummy for feedback before without disrupting their other work.

By using the following command, we can make a new branch so that Dracula and Mummy won't need to worry with Wolfman's changes, until he's ready to show them.

<div class="in" markdown="1">
~~~
git branch shelter_plan
~~~
</div>

We just created a new branch, now we need to switch to it to start working. We do that by doing the following:

<div class="in" markdown="1">
~~~
git checkout shelter_plan
~~~
</div>

If we type `git branch` again, we should see what branch we're on. Let's give that a try. We should get the following:

<div class="out" markdown="1">
~~~
* shelter_plan
  master
~~~
 </div>

Now, it's time to do our git workflow to add the `shelter.txt` file we made earlier to our new branch. This is how we do that:

<div class="in" markdown="1">
~~~
git add shelter.txt
git commit 
git push origin shelter_plan
~~~
</div>

So, we just pushed `shelter.txt` to our new branch of the main project, shelter_plan! Let's see if that worked.

<div class="in" markdown="1">
~~~
git branch
~~~
</div>

We see that there are two branches possible:

<div class="out" markdown="1">
~~~
* shelter_plan
  master
~~~
</div>

We see that we are currently in our "shelter_plan" branch, currently, by observing the `*` next to it. 

The utility of creating a branch serves a few purposes here: 
1. This is a benefit for both Wolfman, as an individual contributor and for Dracula and Mummy as collaborators. 
2. Wolfman has made changes without calling Dracula and Mummy's attention. (Which is good because Dracula and Mummy are rough around the edges!)

Now that our first draft for `shelter.txt` is complete, let's run our floor plan by Dracula and Mummy and get their feeback.

*Introduce Pull Request Here!* 

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
* For Dracula: Coffin 80x80
* For Mummy: Sarcaphogus 80x80 
* For me: Bedroom with no windows 80x80
~~~
</div>

Ah, that's much more fair. Wolfman is now ready to commit these changes to the main project. So he needs to switch branches to do so by using the command `git checkout`.

Let's try that on for size:

<div class="in" markdown="1">
~~~
git checkout master
git branch
~~~
</div>

We should see that we are back to the main project, known as the master in our output, as shown below:

<div class="out" markdown="1">
~~~
* master
  shelter_plan
~~~
</div>

Now that we are in the master branch, we can merge the work we did in our shelter branch to the master, so Dracula and Mummy can suggest changes.

<div class="out" markdown="1">
~~~
git merge shelter_plan
~~~
</div>

Great! Now we can "push our changes to the master branch!"

<div class="in" markdown="1">
~~~
git add shelter.txt
git commit sheter.txt
git push
~~~
</div>

Remember, you can always use `git status` to check in if there's anything in your working directory to commit. If your working directory is clean, git will tell you. 

Let's checkout our repository online to ensure we have added the right dimensions to shelter.txt!
Fill in your github username in the following link to YOUR planets repo as so: https://github.com/yourusername/planets and navigate to shelter.txt. 
