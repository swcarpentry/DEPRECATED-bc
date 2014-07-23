---
layout: lesson
root: ../..
title: Working Remotely
---
Let's take a closer look at what happens when we use the shell
on a desktop or laptop computer.
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
We call this a [remote login](../../gloss.html#remote-login),
and the other computer a remote computer.
Once we do this,
everything we type is passed to a shell running on the remote computer.
That shell runs those commands on our behalf,
just as a local shell would,
then sends back output for our computer to display.

The tool we use to log in remotely is the [secure shell)(../../gloss.html#secure-shell),
or SSH.
In particular, the command `ssh username@computer`
runs SSH and connects to the remote computer we have specified.
After we log in,
we can use the remote shell to use the remote computer's files and directories.
Typing `exit` or Control-D
terminates the remote shell and returns us to our previous shell.

In the example below,
the remote machine's command prompt is `moon>`
instead of just `$`.
To make it clearer which machine is doing what,
we'll indent the commands sent to the remote machine
and their output.

~~~
$ pwd
~~~
{:class="in"}
~~~
/users/vlad
~~~
{:class="out"}
~~~
$ ssh vlad@moon.euphoric.edu
Password: ********
~~~
{:class="in"}
~~~
    moon> hostname
~~~
{:class="in"}
~~~
    moon
~~~
{:class="out"}
~~~
    moon> pwd
~~~
{:class="in"}
~~~
    /home/vlad
~~~
{:class="out"}
~~~
    moon> ls -F
~~~
{:class="in"}
~~~
    bin/     cheese.txt   dark_side/   rocks.cfg
~~~
{:class="out"}
~~~
    moon> exit
~~~
{:class="in"}
~~~
$ pwd
~~~
{:class="in"}
~~~
/users/vlad
~~~
{:class="out"}

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
The command's name combines `cp`'s and `ssh`'s,
and so does its operation.
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
~~~
{:class="in"}
~~~
results.dat              100%  9  1.0 MB/s 00:00
~~~
{:class="out"}

Copying a whole directory is similar:
we just use the `-r` option to signal that we want copying to be recursive.
For example,
this command copies all of our results from the backup server to our laptop:

~~~
$ scp -r vlad@backupserver:backups ./backups
Password: ********
~~~
{:class="in"}
~~~
results-2011-09-18.dat              100%  7  1.0 MB/s 00:00
results-2011-10-04.dat              100%  9  1.0 MB/s 00:00
results-2011-10-28.dat              100%  8  1.0 MB/s 00:00
results-2011-11-11.dat              100%  9  1.0 MB/s 00:00
~~~
{:class="out"}

Here's one more thing SSH can do for us.
Suppose we want to check whether we have already created the file
`backups/results-2011-11-12.dat` on the backup server.
Instead of logging in and then typing `ls`,
we could do this:

~~~
$ ssh vlad@backupserver "ls results"
Password: ********
~~~
{:class="in"}
~~~
results-2011-09-18.dat  results-2011-10-28.dat
results-2011-10-04.dat  results-2011-11-11.dat
~~~
{:class="out"}

SSH takes the argument after our remote username
and passes them to the shell on the remote computer.
(We have to put quotes around it to make it look like a single argument.)
Since those arguments are a legal command,
the remote shell runs `ls results` for us
and sends the output back to our local shell for display.

> #### All Those Passwords
>
> Typing our password over and over again is annoying,
> especially if the commands we want to run remotely are in a loop.
> To remove the need to do this,
> we can create an [SSH key](../../gloss.html#ssh-key)
> to tell the remote machine
> that it should always trust us.
> We discuss SSH keys in our intermediate lessons.
