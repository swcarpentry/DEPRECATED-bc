---
layout: lesson
root: ../..
title: Basic Relational Structure of Git Repositories
level: intermediate
---

##The Bones of Git

In earlier lessons, we've explored some of the basic procedures for creating and navigating a git repository; in this lesson, we'll start learning about the formal structure of a git repo, so that we can understand some of git's more advanced uses when we come to them in later lessons.

As you've already seen, the basic building block of a git repo is a **revision**, each of which has a unique **revision hash** - that long string of random-seeming characters associated with each revision visible when you do a `git log`.  If all git did was create these snapshots of our code, it would still be pretty helpful - but git's true power comes from how it organizes these revisions, in a relational structure.

Git's relational structure is essentially one of inheritance - each revision is a child of at least one parent whose content is merged into the child revision (except of course the very first revision; it's created as the starting point of this little family tree).  For now, let's start thinking about the simplest example of this structure: a simple chain of revisions, where each revision only inherits from the revision that came immediately before it, and is only inherited by the next revision after it:

```
[Rev 0]-----[Rev 1]-----[Rev 2]
```

In earlier lessons, this is exactly the structure you were building; every time you did `git commit`, another revision was getting tacked onto the chain.  In git parlance, this chain is called a **branch**; in later lessons, as the name suggests, we'll explore how these branches can spit up and merge together, like the branches of a plant.  But for now, we're just going to deal with one simple, linear chain of commits.  For the curious, in a new git repo, try doing `git branch`; it should return `* master`, which indicates you are in a branch called `master` - that's the name git creates by default for our first chain of commits.  This is just a curiosity for now, but it'll be useful to understand later.

Git has more structure than just the chain of commits, too.  Every project has something called the `HEAD`, which is a term you've already encountered.  In our simple single linear branch, `HEAD` is a pointer to whatever revision git thinks is supposed to be the parent of the revision that's going to come next.  By default, it always points at the most recent revision in the chain:

```
[Rev 0]-----[Rev 1]-----[Rev 2]
                           ^
                           |
               HEAD---------
```

You've already seen the first simple use of git's relational structure by way of `HEAD` - in the above diagram, asking for `HEAD~1` refers to `[Rev 1]`, and `HEAD~2` refers to `[Rev 0]`; we can simply count backwards from `HEAD` in order to ask for earlier revisions.

As we saw in earlier lessons, just firing off a `git commit` at any old time won't actually do anything if you haven't done a `git add` on anything first.  Every time we do `git add`, we're taking the current state of a file, and adding it to something called the **index**.  The index is like a staging area where git keeps track of the files to be added to the next revision, and we can think of it as living at the end of our chain of revisions, sticking off of whatever revision `HEAD` is pointing at:

```
[Rev 0]-----[Rev 1]-----[Rev 2]-----(Index)
                           ^
                           |
               HEAD---------
``` 

The final piece of structure to git is called the **working tree**, which is actually just a fancy name for whatever you've saved in the directory containing your project.  Every time you change a file and save it, you've modified the working tree.  We can think of it as a box of stuff that lives alongside your git project:

```
[Rev 0]-----[Rev 1]-----[Rev 2]-----(Index)        ----------------------
                           ^                       |    working tree    |
                           |                       ----------------------
               HEAD---------
``` 

Every time you do `git add myFile.txt`, you're taking the current state of `myFile.txt` and adding it to the index, ready to be part of the next revision; every time you do `git commit -m 'my message'`, you're packaging up the contents of the index into the next revision in the chain, and moving `HEAD` to point at that new revision.

Let's step through an imaginary work cycle like you've done before, but now looking at exactly what's going on inside git at each step.  Say you started with just one file in your project, `main.txt`; I'll list it inside the revision boxes like so:

```
[Rev 0: main.txt]-----[Rev 1: main.txt]-----[Rev 2: main.txt]-----(Index)        ----------------------
                                                         ^                       |    working tree    |
                                                         |                       |    main.txt        |
                                                         |                       ----------------------
                                             HEAD---------
``` 

At each revision, you did some changes to `main.txt` and committed them, as usual.  Note some version of `main.txt` is hanging out in the working tree too - it may or may not be different from its earlier versions.  Now let's say you create a new file, `new.txt`, and save it in your project directory.  That means it's just in the working tree, and not in your repo yet:

```
[Rev 0: main.txt]-----[Rev 1: main.txt]-----[Rev 2: main.txt]-----(Index)        ----------------------
                                                         ^                       |    working tree    |
                                                         |                       |    main.txt        |
                                                         |                       |    new.txt         |
                                                         |                       ----------------------
                                             HEAD---------
``` 

Now let's add it to the index with `git add new.txt`:

```
[Rev 0: main.txt]-----[Rev 1: main.txt]-----[Rev 2: main.txt]-----(Index: new.txt)        ----------------------
                                                         ^                                |    working tree    |
                                                         |                                |    main.txt        |
                                                         |                                |    new.txt         |
                                                         |                                ----------------------
                                             HEAD---------
```

And finally, let's commit that as a new revision, with `git commit -m 'a great new work'`

```
[Rev 0: main.txt]-----[Rev 1: main.txt]-----[Rev 2: main.txt]-----[Rev 3: main.txt -----(Index)
                                                                          new.txt]                     ----------------------
                                                                            ^                          |    working tree    |
                                                                            |                          |    main.txt        |
                                                                            |                          |    new.txt         |
                                                                            |                          ----------------------
                                                                HEAD---------
```

A new revision has been created which records the changes you staged to new.txt (in this case, 'changes' being the creation of the file), the index has been emptied, and `HEAD` now points at the latest revision, `[Rev 3]`.  That's all there is to it!  That's what git is doing under the hood when you use the basic tools you've learned up to this point.

Now that we understand this basic relational structure, we're ready to learn about a new tool for walking around our chain of commits: `git reset`.

##`git reset`

Let's suppose that in the course of our work, we make some kind of mistake, and we want to rewind the project to an earlier point - `git reset` might be the right tool for the job; think of it as git's big 'undo' button.  `git reset` has three options to choose from: `--hard`, `--mixed`, or `--soft`.  The simplest to understand is the most drastic: `git reset --hard HEAD` will empty the index, change everything in the working tree to match the version recorded at `HEAD`, and move the `HEAD` pointer to point at whatever revision we pointed at (which in the example of `git reset --hard HEAD` of course doesn't move `HEAD` at all, but it would if we did something like `git reset --hard HEAD~1`, or used an earlier revision hash instead of `HEAD`).  So continuing our example from the last section, if we did `git reset --hard HEAD~1`, our repo would look something like:

```
[Rev 0: main.txt]-----[Rev 1: main.txt]-----[Rev 2: main.txt]-----[Rev 3: main.txt -----(Index)
                                                    ^                     new.txt]                     ----------------------
                                                    |                                                  |    working tree    |
                                                    |                                                  |    main.txt        |
                                                    |                                                  ______________________
                                                    |                                                 
                                        HEAD---------
```

Notice something very crucial about `git reset --hard`: `new.txt` has vanished from the working tree (and, though it's harder to depict, any changes made to `main.txt` in `[Rev 3]` are undone in the working tree, too); this means that **reset --hard can destroy work if you're not careful!**  Luckily for us, we committed `new.txt` in `[Rev 3]` before we reset, and `[Rev 3]` still exists; we could do `git reset --hard <Rev 3's hash>` to go back to `[Rev 3]`; `reset` can go forward as well as backwards, if we have the hash around to tell it where to go.  Finally, note that `HEAD` is now pointing at the revision we reset to; this means that if we start committing again from this point, our chain is going to carry on from `[Rev 2]`, with `[Rev 3]` now out of the sequence traversed by `HEAD~1`, `HEAD~2`... : 

```
[Rev 0: main.txt]-----[Rev 1: main.txt]-----[Rev 2: main.txt]-----[Rev 3: main.txt
                                                    |                     new.txt]                     ----------------------
                                                    |                                                  |    working tree    |
                                                    |                                                  |    main.txt        |
                                                    |-----[Rev 4: main.txt]-----(Index)                ______________________
                                                                 ^
                                                                 |                                                 
                                                     HEAD---------
```

As promised, there are less drastic things you can do than `git reset --hard`.  Next up is `git reset --mixed <revision>`, which empties the index and moves `HEAD` to point at `<revision>` just like `reset --hard`, but doesn't touch the working tree - that way, your working tree is totally safe and no work living there will be destroyed.  A common use for this would be `git reset --mixed HEAD`, which doesn't move the `HEAD` pointer, but just empties the index so you can commit things differently.

The final and most tame option for reset is `git reset --soft <revision>`.  The `--soft` flag tells reset to leave the working tree *and* the index alone, and just move the `HEAD` pointer to `<revision>`, so we can use an earlier commit as the parent for our next one.  Notice that `git reset --soft HEAD` actually does nothing - this one only makes sense for revisions other than `HEAD`.  A common use for the `--soft` flag is if you've just made a commit, and realize you want to change the commit message; just go `git reset --soft HEAD~1`, and then recommit with a new message.

To summarize all the flavours of commit and their popular uses, here's a table:

| Command             | Effect                                         | Common Use                                 |
|---------------------|------------------------------------------------|--------------------------------------------|
| `git reset --soft`  | No changes to working tree or index            |Change the commit message you just made     |
| `git reset --mixed` | Remove staged changes from index               |Change what's going into the next commit    |
| `git reset --hard`  | Remove changed files in index and working tree |Completely abandon everything and go back   |
