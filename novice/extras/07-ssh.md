---
layout: lesson
root: ../..
title: Working Remotely
---
<div class="objectives" markdown="1">

#### Objectives
*   Explain what is SSH
*   Explain what an SSH key is
*   Generate your own SSH key pair
*   Add your SSH key to an remote server
*   Learn how to use your SSH key

</div>

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

The tool we use to log in remotely is the [secure shell](../../gloss.html#secure-shell),
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

### SSH Keys

Typing our password over and over again is annoying,
especially if the commands we want to run remotely are in a loop.
To remove the need to do this,
we can create an [SSH key](../../gloss.html#ssh-key)
to tell the remote machine
that it should always trust us.

SSH keys come in pairs, a public key that gets shared with services like GitHub,
and a private key that is stored only on your computer. If the keys match,
you're granted access.

The cryptography behind SSH keys ensures that no one can reverse engineer your
private key from the public one.

The first step in using SSH authorization is to generate your own key pair.

You might already have an SSH key pair on your machine. You can check to see if
one exists by moving to your `.ssh` directory and listing the contents.

~~~
$ cd ~/.ssh
$ ls
~~~
{:class="in"}

If you see `id_rsa.pub`, you already have a key pair and don't need to create a
new one.

If you don't see `id_rsa.pub`, use the following command to generate a new key
pair. Make sure to replace `your@email.com` with your own email address.

~~~
$ ssh-keygen -t rsa -C "your@email.com"
~~~
{:class="in"}

When asked where to save the new key, hit enter to accept the default location.

~~~
Generating public/private rsa key pair.
Enter file in which to save the key (/Users/username/.ssh/id_rsa):
~~~
{:class="out"}

You will then be asked to provide an optional passphrase. This can be used to
make your key even more secure, but if what you want is avoiding type your
password every time you can skip it by hitting enter twice.

~~~
Enter passphrase (empty for no passphrase):
Enter same passphrase again:
~~~
{:class="out"}

When the key generation is complete, you should see the following confirmation:

~~~
Your identification has been saved in /Users/username/.ssh/id_rsa.
Your public key has been saved in /Users/username/.ssh/id_rsa.pub.
The key fingerprint is:
01:0f:f4:3b:ca:85:d6:17:a1:7d:f0:68:9d:f0:a2:db your@email.com
The key's randomart image is:
+--[ RSA 2048]----+
|                 |
|                 |
|        . E +    |
|       . o = .   |
|      . S =   o  |
|       o.O . o   |
|       o .+ .    |
|      . o+..     |
|       .+=o      |
+-----------------+
~~~
{:class="out"}

The random art image is an alternate way to match keys but we won't be needing this.

Now you need to send your public key to the server you want to connect. Display
the contents of your new public key file with `cat`:

~~~
$ cat ~/.ssh/id_rsa.pub
~~~
{:class="in"}
~~~
ssh-rsa AAAAB3NzaC1yc2EAAAABIwAAAQEA879BJGYlPTLIuc9/R5MYiN4yc/YiCLcdBpSdzgK9Dt0Bkfe3rSz5cPm4wmehdE7GkVFXrBJ2YHqPLuM1yx1AUxIebpwlIl9f/aUHOts9eVnVh4NztPy0iSU/Sv0b2ODQQvcy2vYcujlorscl8JjAgfWsO3W4iGEe6QwBpVomcME8IU35v5VbylM9ORQa6wvZMVrPECBvwItTY8cPWH3MGZiK/74eHbSLKA4PY3gM4GHI450Nie16yggEg2aTQfWA1rry9JYWEoHS9pJ1dnLqZU3k/8OWgqJrilwSoC5rGjgp93iu0H8T6+mEHGRQe84Nk1y5lESSWIbn6P636Bl3uQ== your@email.com
~~~
{:class="out"}

Copy the contents of the output. Login to the server you want to connect using
your SSH keys.

~~~
$ ssh vlad@moon.euphoric.edu
Password: ********
~~~
{:class="in"}

Paste the content that you copy at the end of `~/.ssh/authorized_keys`.

~~~
    moon> nano ~/.ssh/authorized_keys`.
~~~
{:class="in"}

After append the content, logout of the remote machine and try login again. If
you setup your SSH key correctly you won't need to type your password.

~~~
    moon> exit
~~~
{:class="in"}
~~~
$ ssh vlad@moon.euphoric.edu
~~~
{:class="in"}

<div class="keypoints" markdown="1">

#### Key Points
*  SSH is a secure alternative to username/password authorization
*  SSH keys are generated in public/private pairs. Your public key can be shared
   with others. The private keys stays on your machine only.
</div>
