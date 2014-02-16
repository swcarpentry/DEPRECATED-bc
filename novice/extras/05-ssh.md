---
layout: lesson
root: ../..
title: Working Remotely
level: novice
---
Let's take a closer look at what happens when we use a desktop or laptop computer.
The first step is to log in
so that the operating system knows who we are and what we're allowed to do.
We do this by typing our username and password;
the operating system checks those values against its records,
and if they match,
runs a shell for us.

As we type commands,
the 1's and 0's that represent the characters we're typing are sent from the keyboard to the shell.
The shell displays those characters on the screen to represent what we type,
and then,
if what we typed was a command,
the shell executes it and displays its output (if any).

What if we want to run some commands on another machine,
such as the server in the basement that manages our database of experimental results?
To do this,
we have to first log in to that machine.
We call this a [remote login](../gloss.html#remote-login),
and the other computer a remote computer.
Once we do this,
everything we type is passed to a shell running on the remote computer.
That shell interacts runs those commands on our behalf,
just as a local shell would,
then sends back output for our computer to display.

The tool we use to log in remotely is the [secure shell)(../gloss.html#secure-shell),
or SSH.
In particular, the command `ssh username@computer`
runs SSH and connects to the remote computer we have specified.
After we log in,
we can use the remote shell to use the remote computer's files and directories.
Typing `exit` or Control-D
terminates the remote shell and returns us to our previous shell.
In the example below, we use highlighting to show our interaction with the remote shell.
We can also see that the remote machine's command prompt is `moon>`
instead of just `$`:

~~~
$ pwd
/users/vlad
$ ssh vlad@moon.euphoric.edu
Password: ********
moon> hostname
moon
moon> pwd
/home/vlad
moon> ls -F
bin/     cheese.txt   dark_side/   rocks.cfg
moon> exit
$ pwd
/users/vlad
~~~

The secure shell is called "secure" to contrast it with an older program called `rsh`,
which stood for "remote shell".
Back in the day,
when everyone trusted each other and knew every chip in their computer by its first name,
people didn't encrypt anything except the most sensitive information when sending it over a network.
However,
that meant that villains could watch network traffic,
steal usernames and passwords,
and use them for all manner of nefarious purposes.
SSH was invented to prevent this (or at least slow it down).
It uses several sophisticated, and heavily tested, encryption protocols
to ensure that outsiders can't see what's in the messages
going back and forth between different computers.

`ssh` has a companion program called `scp`,
which stands for "secure copy".
It allows us to copy files to or from a remote computer using the same kind of connection as SSH.
The syntax is a simple mix of `cp`'s and `ssh`'s.
To copy a file,
we specify the source and destination paths,
either of which may include computer names.
If we leave out a computer name,
`scp` assumes we mean the machine we're running on.
For example,
this command copies our latest results to the backup server in the basement,
printing out its progress as it does so:

~~~
$ scp results.dat vlad@backupserver:backups/results-2011-11-11.dat
Password: ********
results.dat              100%  9  1.0 MB/s 00:00
~~~

Copying a whole directory is similar:
we just use the `-r` option to signal that we want copying to be recursive.
For example,
this command copies all of our results from the backup server to our laptop:

~~~
$ scp -r vlad@backupserver:backups ./backups
Password: ********
results-2011-09-18.dat              100%  7  1.0 MB/s 00:00
results-2011-10-04.dat              100%  9  1.0 MB/s 00:00
results-2011-10-28.dat              100%  8  1.0 MB/s 00:00
results-2011-11-11.dat              100%  9  1.0 MB/s 00:00
~~~

To close off this chapter,
let's look at one more thing SSH can do for us.
Suppose we want to check whether we have already created the file
`backups/results-2011-11-12.dat` on the backup server.
Instead of logging in and then typing `ls`,
we could do this:

~~~
$ ssh vlad@backupserver ls results
Password: ********
results-2011-09-18.dat  results-2011-10-28.dat
results-2011-10-04.dat  results-2011-11-11.dat
~~~

SSH has taken the arguments after our username and the name of the computer we want to run on
and passed them to the shell on the remote computer.
Since those arguments are a legal command,
the remote shell has run `ls results` for us
and sent the output back to our local shell for display.
