## Getting Started.

In these exercises you will get to practice using a version control system,
Subversion.  In these exercises we will be using the SmartSVN client because it
is cross-platform.  Download and install the SmartSVN GUI from
[here](http://www.syntevo.com/smartsvn/index.html).  You may also choose to use
a [command-line client](http://subversion.apache.org/), or another Subversion
client.  If you do so, you'll have to adapt the solutions we give to your client
since our solutions are for the SmartSVN client.

Please attempt these exercises by looking back at what you learned in the
Version Control lecture.  If you get stuck, or don't know what to do, there are
links to solutions at the bottom of each exercise.


## Exercise 1: Checking out.

We've set up a repository for the class.  The repository URL is:

    http://svn.software-carpentry.org/spring2011

You will need a username and password to access the repository.  Your _username_
is your family name followed by the first letter of your first name, all lower
case with spaces removed.  So, "John von Neumann" would have the username
"`vonneumannj`".  Your _password_ is "swc2011".

Use SmartSVN to check out a working copy of your repository.

**Answer**:
> 1. From the "Project" menu in the upper left corner, select "Check out".
> 
> 2. Select the radio button for "Quick Checkout".  Fill in the "URL" field with
> the repository URL we gave above.
> 
> 3. On the "Local Directory" line, click the "Browse" button, and choose a
> place where you want to store the working copy of the repository on your
> computer.
> 
> 4. Click "Next".
> 
> 5. If you get asked "Do you really want to check out the whole repository?",
> select "Check out".
> 
> 6. Click "Finish".


## Exercise 2: Modifying and Committing.

In your working copy you should see a folder called "participants" containing a
folder with the username of each member of the class, including your own. You
should find two files in your folder: `planets.txt` and `temperature.txt`.

1. Open `temperature.txt`, and change the "January Low" for Winnipeg to -23.5.

2. Save the file.

3. Commit your change to the repository, being sure to include an informative
comment.

4. When prompted for authentication, your user name is, again, your last name,
followed by the first letter of your first name.  Your password is _swc2011_.

**Answer:**	
> 1. Double click `temperature.txt` to open it in your text editor.
> 
> 2. Change the "January Low" for Winnipeg to -23.5.
> 
> 3. Save the file
> 
> 4. Back in SmartSVN, select the `temperature.txt` file, and click "Commmit"
> 
> 5. In the "Commit Message" box, give a meaningful comment, like "Changed the
> January Low for Winnipeg to -23.5."
> 
> 6. Click the "Commit" button
> 
> 7. You will be asked for authentication.  Your username is the your last name,
> followed by the first letter of your first name.  Your password is _swc2011_.
> 
> 8. Click "Login".


## Exercise 3: Modifying and Committing... Again.

1. Open `planets.txt`, and add a new second line, with _* 10^6 miles_ typed
directly under _Mean Distance From Sun_.  (See the picture below).

2. Save the file.

3. Commit your change to the repository, being sure to include an informative
comment.

**Answer:**
> We're not helping you this time!  This is the exact same process as Exercise 2.


## Exercise 4: Resolve a Conflict

In this exercise, a different user is going to make a change to your
`planets.txt` file and check it in before you check in your changes to the same
file.  This will mean you will have to resolve the conflicting lines before you
can check in your changes. 

We have written a special program that modifies the `planets.txt` in your
repository folder as if it were another user.  Visit the following URL: 

    http://svn.software-carpentry.org/exercise_4.cgi

and type in your username.  When you hit submit you should receive a message
saying `planets.txt` was successfully modified. DO NOT UPDATE YOUR WORKING COPY
YET.  Please follow the instructions below. 

1. Open up your copy of `planets.txt`
2. There's a typo!  Change "Marz" to "Mars".
3. _* 10^6 miles_ isn't very readable.  Change it to _million miles_.
4. Save the file.
5. Try to commit your changes.  You should get an error.
6. Update, resolve any conflicts, and commit your changes.

**Answer:**
> 1. Open SmartSVN.  Select your working copy from the "Open Existing Project"
> pane.  Click "Okay".
> 2. Double click `planets.txt` to open it in your text editor.
> 3. Make the changes we described, and save the file.
> 4. Back in SmartSVN, select the `planets.txt` file, and click the "Commit" at
> the top. Type your message, and click "Commit".
> 5. You will get an error, since you are not at the most recent revision.  Click
> "Okay".
> [![](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex06_01-300x120.jpg)](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex06_01.jpg)
> 6. Click the "Update" button at the top.
> 7. Now you need to resolve the conflict.  Right-click on `planets.txt`.  Go to
> "Query" and then "Conflict Solver".
> [![](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_06_02-300x113.jpg)](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_06_02.jpg)
> 8. In the Conflict Solver, the middle pane shows the current state of
> `planets.txt`.  We want to select the "million miles" instead of the "1000000
> miles", so click the double arrows (>>) next to the "million miles" line in the
> left pane.  In the middle pane, we can see that our "Mars" correction has been
> merged, and so has the information about Saturn that the other user added.
> These changes are fine, so we just leave them. Click "Save" at the top, and then
> close the Conflict Solver.
> [![](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_06_03-300x162.jpg)](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_06_03.jpg)
> 9. Select `planets.txt`, and click the "Commit" button at the top.  This time the commit should be successful.



## Exercise 5: Reversing Local Changes.
1. Open `temperature.txt`, and change the January high for Toronto to +600.
2. Save the file, and close your text editor.
3. Now use SVN to undo the change you made.

**Solution:**
> 1. Open `temperature.txt`, make the changes, save the file, and close your
> text editor.
> 2. In SmartSVN, right-click on `temperature.txt`, and select "Revert".
> 3. If it asks you if you are sure, choose "Revert".


## Exercise 6: Reversing Committed Changes.
	
1. Again, open `temperature.txt`, and change the January high for Toronto to
+600.

2. Save the file, and close your text editor.

3. Commit your changes.

4. Use SVN to reverse the changes in the latest version of the repository.

**Solution:**
> 1. Open `temperature.txt`, make the changes, save the file, and close your text
> editor.
> 
> 2. In SmartSVN, commit the changes.
> 
> 3. From the "Modify" menu at the top, select "Merge".
> [![](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_08_01-300x195.jpg)](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_08_01.jpg)
> 
> 4. Type in the revision range.  The start value for the range should be the
> revision before your most recent commit.  The end value of the range should be
> HEAD.  If you've been following these exercises to the letter, the range should
> be "6-HEAD".  If you've done some extra commits along the way, your start value
> might be more than 6.
> 
> 5. Select the "Reverse Merge" option.  (You might have to deselect the "Enable
> xMerge" option first.)
> 
> 6. Click the "Merge" button.
> 
> 7. Open `temperature.txt` and verify that the change has been reversed.
> 
> 8. Back in SmartSVN, click the "Commit" button at the top, and in the message,
> explain that you have changed the Toronto January high back to -1.1.  Click
> "Commit".


## Exercise 7: Adding and Committing a New File.

You will now make a new file, and add it to the repository.

1. Make a new file called `interests.txt` in your personal repository folder.
        
2. This file should contain a single line that starts with your name, followed
by a colon, followed by a comma-separated list of your research interests.  For
example:

    Henrietta Leavitt: astronomy, Cepheid variables, period-luminosity
    relationship

3. Add this file to the repository.

4. Commit your change to the repository, being sure to include an informative
comment.

**Answer:**
> 1. Open your favourite text editor, and type your name and research interests
> into in the format given above.
> 
> 2. Save this file as `interests.txt` in the working directory of the
> repository (this is the same directory you used selected for "Local Directory"
> when you first checked out.)
> 
> 3. If you've forgotten where your working copy is, take a look at the
> Directories panel in SmartSVN.
> [![](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_04_01-300x117.jpg)](http://software-carpentry.org/blog/wp-content/uploads/2010/10/ex_04_01.jpg)
> 
> 4. Go back to SmartSVN.  You should see `interests.txt` in the main panel.  If
> not, click "View" in the top menu, and go to "Refresh".
> 
> 5. Select `interests.txt` and click the "Add" button.
> 
> 6. Click the "Commit" button (next to the "Add" button)
> 
> 7. Put a meaningful comment in the "Commit Message" box, and click "Commit."


## Exercise 8: Editing a shared file.

In this exercise, the unexpected may happen.  You are going to add what you
wrote in your personal `interests.txt` file in Exercise 7, to an `interests.txt`
file shared by the rest of the class.  Depending on who else is working on these
exercises at the same time as you are, you might find you need to merge your
changes in with theirs.  Who knows!

1. Open up the `interests.txt` file that is in the top-level folder of the
repository (i.e. at the same level as the `participants` folder.

2. Add your name and interests on a new line, as you did in Exercise 7.

3. Commit your changes with a useful log message. You may need to update and
resolve conflicts before you do.
