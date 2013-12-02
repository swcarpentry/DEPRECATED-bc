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

**FIXME**: <a name="additive-color-model"></a>

**FIXME**: <a name="aggregation-function"></a>

**FIXME**: <a name="alias"></a>

**FIXME**: <a name="assertion"></a>

**FIXME**: <a name="assignment"></a>

**FIXME**: <a name="atomic-value"></a>

**FIXME**: <a name="cascading-delete"></a>

**FIXME**: <a name="case-insensitive"></a>

**FIXME**: <a name="change-set"></a>

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

**FIXME**: <a name="comment"></a>

**FIXME**: <a name="conditional-statement"></a>

**FIXME**: <a name="confirmation-bias"></a>

**FIXME**: <a name="conflict"></a>

**FIXME**: <a name="cross-product"></a>

**FIXME**: <a name="csv"></a>

**current working directory**: <a name="current-working-directory"></a>
The directory that [relative paths](#relative-path) are calculated from;
equivalently,
the place where files referenced by name only are searched for.
Every [process](#process) has a current working directory.
The current working directory is usually referred to using the shorthand notation `.` (pronounced "dot").

**FIXME**: <a name="cursor"></a>

**FIXME**: <a name="database-manager"></a>

**FIXME**: <a name="default-parameter-value"></a>

**FIXME**: <a name="defensive-programming"></a>

**FIXME**: <a name="delimiter"></a>

**disk blocks**: <a name="disk-block"></a>
FIXME

**FIXME**: <a name="docstring"></a>

**FIXME**: <a name="documentation"></a>

**FIXME**: <a name="dotted-notation"></a>

**FIXME**: <a name="empty-string"></a>

**FIXME**: <a name="encapsulation"></a>

**FIXME**: <a name="field-database"></a>

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

**filesystem**: <a name="filesystem"></a>
A set of files, directories, and I/O devices (such as keyboards and screens).
A filesystem may be spread across many physical devices,
or many filesystems may be stored on a single physical device;
the [operating system](#operating-system) manages access.

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

**FIXME**: <a name="float"></a>

**FIXME**: <a name="foreign-key"></a>

**FIXME**: <a name="function-composition"></a>

**graphical user interface** (GUI): <a name="gui"></a>
FIXME

**home directory**: <a name="home-directory"></a>
FIXME

**FIXME**: <a name="immutable"></a>

**FIXME**: <a name="import"></a>

**FIXME**: <a name="in-place-operator"></a>

**FIXME**: <a name="index"></a>

**FIXME**: <a name="infective-license"></a>

**FIXME**: <a name="inner-loop"></a>

**FIXME**: <a name="integer"></a>

**FIXME**: <a name="invariant"></a>

**loop**: <a name="for-loop"></a>
FIXME

**loop body**: <a name="loop-body"></a>
FIXME

**FIXME**: <a name="loop-variable"></a>

**FIXME**: <a name="member"></a>

**FIXME**: <a name="method"></a>

**FIXME**: <a name="mutable"></a>

**FIXME**: <a name="notional-machine"></a>

**orthogonality**: <a name="orthogonality"></a>
FIXME

**FIXME**: <a name="outer-loop"></a>

**FIXME**: <a name="parameter"></a>

**parent directory**: <a name="parent-directory"></a>
The directory that "contains" the one in question.
Every directory in a file system except the [root directory](#root-directory) has a parent.
A directory's parent is usually referred to using the shorthand notation `..` (pronounced "dot dot").

**pipe**: <a name="pipe"></a>
A connection from the output of one program to the input of another.
When two or more programs are connected in this way, they are called a "pipeline".

**pipe and filter**: <a name="pipe-and-filter"></a>
FIXME

**FIXME**: <a name="postcondition"></a>

**FIXME**: <a name="precondition"></a>

**FIXME**: <a name="prepared-statement"></a>

**FIXME**: <a name="primary-key"></a>

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

**FIXME**: <a name="query"></a>

**FIXME**: <a name="range"></a>

**read-evaluate-print loop** (REPL): <a name="repl"></a>
FIXME

**FIXME**: <a name="record-database"></a>

**redirection**: <a name="redirection"></a>
FIXME

**FIXME**: <a name="referential-integrity"></a>

**FIXME**: <a name="regression"></a>

**regular expressions** (RE): <a name="regular-expression"></a>
A pattern that specifies a set of character strings.
REs are most often used to find sequences of characters in strings.

**FIXME**: <a name="relational-database"></a>

**relative path**: <a name="relative-path"></a>
A [path](#path) that specifies the location of a file or directory
with respect to the [current working directory](#current-working-directory).
Any path that does not begin with a separator character ("/" or "\") is a relative path.
See also: [absolute path](#absolute-path).

**FIXME**: <a name="repository"></a>

**FIXME**: <a name="repository-clone"></a>

**FIXME**: <a name="repository-merge"></a>

**FIXME**: <a name="repository-remote"></a>

**FIXME**: <a name="resolve"></a>

**FIXME**: <a name="return-statement"></a>

**FIXME**: <a name="rgb"></a>

**root directory**: <a name="root-directory"></a>
The top-most directory in a [filesystem](#filesystem).
Its name is "/" on Unix (including Linux and Mac OS X) and "\" on Microsoft Windows.

**FIXME**: <a name="sentinel-value"></a>

**FIXME**: <a name="shape"></a>

**shell quoting**: <a name="shell-quoting"></a>
FIXME

**shell script**: <a name="shell-script"></a>
FIXME

**FIXME**: <a name="silent-failure"></a>

**FIXME**: <a name="slice"></a>

**FIXME**: <a name="sql"></a>

**FIXME**: <a name="sql-injection"></a>

**FIXME**: <a name="stack-frame"></a>

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

**FIXME**: <a name="stride"></a>

**FIXME**: <a name="string"></a>

**sub-directories**: <a name="sub-directory"></a>
FIXME

**tab completion**: <a name="tab-completion"></a>
FIXME

**FIXME**: <a name="table"></a>

**FIXME**: <a name="test-driven-development"></a>

**FIXME**: <a name="test-oracle"></a>

**FIXME**: <a name="tuple"></a>

**FIXME**: <a name="type"></a>

**variable**: <a name="variable"></a>
FIXME

**FIXME**: <a name="version-control"></a>

**wildcard**: <a name="wildcard"></a>
A character used in pattern matching.
In the Unix shell,
the wildcard "*" matches zero or more characters,
so that `*.txt` matches all files whose names end in `.txt`.
