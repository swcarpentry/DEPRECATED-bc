---
---
Glossary
========

**absolute path**: <a name="absolute-path"></a>
A [path](#path) that refers to a particular location in a file system.
Absolute paths are usually written with respect to the file system's
[root directory](#root-directory),
and begin with either "/" (on Unix) or "\" (on Microsoft Windows).
See also: [relative path](#relative-path).

**command shell**: <a name="shell"></a>
A command-line user interface program, such as Bash (the Bourne-Again Shell) or the Microsoft Windows DOS shell.
Shells commonly execute a [read-evaluate-print loop](#repl):
when the user enters a command in response to a prompt,
the shell either executes the command itself,
or runs the program that the command has specified.
In either case,
output is sent to the shell window
and the user is prompted to enter another command.

**command-line interface** (CLI): <a name="cli"></a>
FIXME

**current working directory**: <a name="current-working-directory"></a>
The directory that [relative paths](#relative-path) are calculated from;
equivalently,
the place where files referenced by name only are searched for.
Every [process](#process) has a current working directory.
The current working directory is usually referred to using the shorthand notation `.` (pronounced "dot").

**disk blocks**: <a name="disk-block"></a>
FIXME

**filesystem**: <a name="filesystem"></a>
A set of files, directories, and I/O devices (such as keyboards and screens).
A filesystem may be spread across many physical devices,
or many filesystems may be stored on a single physical device;
the [operating system](#operating-system) manages access.

**filename extension**: <a name="filename-extension"></a>
The portion of a file's name that comes after the final "." character.
By convention this identifies the file's type:
`.txt` means "text file",
`.png` means "Portable Network Graphics file",
and so on.
These conventions are not enforced by most operating systems:
it is perfectly possible to name an MP3 sound file `homepage.html`.
Since many applications use filename extensions to identify the [MIME type](#mime-type) of the file,
misnaming files may cause those applications to fail.

**filter**: <a name="filter"></a>
A program that transforms a stream of data.
Many Unix command-line tools are written as filters:
they read data from [standard input](#standard-input),
process it,
and write the result to [standard output](#standard-output).

**flag**: <a name="command-line-flag"></a>
A terse way to specify an option or setting to a command-line program.
By convention Unix applications use a dash followed by a single letter,
such as `-v`,
or two dashes followed by a word,
such as `--verbose`,
while DOS applications use a slash,
such as `/V`.
Depending on the application, a flag may be followed by a single argument, as in `-o /tmp/output.txt`.

**graphical user interface** (GUI): <a name="gui"></a>
FIXME

**home directory**: <a name="home-directory"></a>
FIXME

**loop body**: <a name="loop-body"></a>
FIXME

**loop**: <a name="for-loop"></a>
FIXME

**orthogonality**: <a name="orthogonality"></a>
FIXME

**parent directory**: <a name="parent-directory"></a>
The directory that "contains" the one in question.
Every directory in a file system except the [root directory](#root-directory) has a parent.
A directory's parent is usually referred to using the shorthand notation `..` (pronounced "dot dot").

**pipe**: <a name="pipe"></a>
A connection from the output of one program to the input of another.
When two or more programs are connected in this way, they are called a "pipeline".

**pipe and filter**: <a name="pipe-and-filter"></a>
FIXME

**process**: <a name="process"></a>
A running instance of a program,
containing code,
variable values,
open files and network connections,
and so on.
Processes are the "actors" that the [operating system](#operating-system) manages;
it typically runs each process for a few milliseconds at a time
to give the impression that they are executing simultaneously.

**prompt**: <a name="prompt"></a>
FIXME

**read-evaluate-print loop** (REPL): <a name="repl"></a>
FIXME

**redirection**: <a name="redirection"></a>
FIXME

**regular expressions** (RE): <a name="regular-expression"></a>
A pattern that specifies a set of character strings.
REs are most often used to find sequences of characters in strings.

**relative path**: <a name="relative-path"></a>
A [path](#path) that specifies the location of a file or directory
with respect to the [current working directory](#current-working-directory).
Any path that does not begin with a separator character ("/" or "\") is a relative path.
See also: [absolute path](#absolute-path).

**root directory**: <a name="root-directory"></a>
The top-most directory in a [filesystem](#filesystem).
Its name is "/" on Unix (including Linux and Mac OS X) and "\" on Microsoft Windows.

**shell script**: <a name="shell-script"></a>
FIXME

**standard input** (stdin): <a name="standard-input"></a>
A process's default input stream.
In interactive command-line applications,
it is typically connected to the keyboard;;
in a [pipe](#pipe),
it receives data from the [standard output](#standard-output) of the preceding process.

**standard output** (stdout): <a name="standard-output"></a>
A process's default output stream.
In interactive command-line applications,
data sent to standard output is displayed on the screen;
in a [pipe](#pipe),
it is passed to the [standard input](#standard-input) of the next process.

**sub-directories**: <a name="sub-directory"></a>
FIXME

**tab completion**: <a name="tab-completion"></a>
FIXME

**variable**: <a name="variable"></a>
FIXME

**wildcard**: <a name="wildcard"></a>
A character used in pattern matching.
In the Unix shell,
the wildcard "*" matches zero or more characters,
so that `*.txt` matches all files whose names end in `.txt`.
