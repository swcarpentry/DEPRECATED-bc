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

**aggregation function**: <a name="aggregation-function"></a>
A function such as `sum` or `max` that combines many values to produce a single result.

**FIXME**: <a name="alias-library"></a>

**assertion**: <a name="assertion"></a>
An expression which is supposed to be true at a particular point in a program.
Programmers typically put assertions in their code to check for errors;
if the assertion fails (i.e., if the expression evaluates as false),
the program halts and produces an error message.
See also: [invariant](#invariant), [precondition](#precondition), [postcondition](#postcondition),

**FIXME**: <a name="assignment"></a>

**atomic value**: <a name="atomic-value"></a>
A value that cannot be decomposed into smaller pieces.
For example,
the number 12 is usually considered atomic
(unless we are teaching addition to school children,
in which case we might decompose it into tens and ones).

**FIXME**: <a name="call-stack"></a>

**cascading delete**: <a name="cascading-delete"></a>
The practice of automatically deleting things in a database
that depend on a record
when that record is deleted.
See also: [referential integrity](#referential-integrity).

**case insensitive**: <a name="case-insensitive"></a>
Treating text as if upper and lower case characters were the same.
See also: [case sensitive](#case-sensitive).

**change set**: <a name="change-set"></a>
A group of changes to one or more files
that are [committed](#commit) to a [version control](#version-control) [repository](#repository)
in a single operation.

**comma-separated values** (CSV): <a name="csv"></a>
A common textual representation for tables
in which the values in each row are separated by commas.

**command-line interface** (CLI): <a name="cli"></a>
An interface based on typing commands,
usually at a [REPL](#repl).
See also: (graphical user interface)[#gui].

**FIXME**: <a name="comment"></a>

**FIXME**: <a name="conditional-statement"></a>

**FIXME**: <a name="confirmation-bias"></a>

**conflict**: <a name="conflict"></a>
A change made by one user of a (version control system)[#version-control]
that is incompatible with changes made by other users.
Helping users (resolve)[#resolve] conflicts
is one of version control's major tasks.

**cross product**: <a name="cross-product"></a>
A pairing of all elements of one set with all elements of another.

**current working directory**: <a name="current-working-directory"></a>
The directory that [relative paths](#relative-path) are calculated from;
equivalently,
the place where files referenced by name only are searched for.
Every [process](#process) has a current working directory.
The current working directory is usually referred to using the shorthand notation `.` (pronounced "dot").

**cursor**: <a name="cursor"></a>
A pointer into a database that keeps track of outstanding operations.

**database manager**: <a name="database-manager"></a>
A program that manages a (relational database)[#relational-database].

**default parameter value**: <a name="default-parameter-value"></a>
A value to use for a parameter if nothing is specified explicitly.

**defensive programming**: <a name="defensive-programming"></a>
The practice of writing programs that check their own operation to catch errors as early as possible.

**delimiter**: <a name="delimiter"></a>
A character or characters used to separate individual values,
such as the commas between columns in a [CSV](#csv) file.

**FIXME**: <a name="disk-block"></a>

**docstring**: <a name="docstring"></a>
Short for "documentation string",
this refers to textual documentation embedded in Python programs.
Unlike comments,
docstrings are preserved in the running program
and can be examined in interactive sessions.

**FIXME**: <a name="documentation"></a>

**FIXME**: <a name="dotted-notation"></a>

**empty string**: <a name="empty-string"></a>
A character string containing no characters,
often thought of as the "zero" of text.

**encapsulation**: <a name="encapsulation"></a>
The practice of hiding something's implementation details
so that the rest of a program can worry about *what* it does
rather than *how* it does it.

**field** (of a database): <a name="field-database"></a>
A set of data values of a particular type,
one for each [record](#record-database) in a (table)[#table-database].

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

**floating point number** (float): <a name="float"></a>
A number containing a fractional part and an exponent.
See also: [integer](#integer).

**for loop**: <a name="for-loop"></a>
A loop that is executed once for each value in some kind of set, list, or range.
See also: [while loop](#while-loop).

**foreign key**: <a name="foreign-key"></a>
One or more values in a (database table)[#table-database]
that identify a (records)[#record-database] in another table.

**FIXME**: <a name="function-composition"></a>

**graphical user interface** (GUI): <a name="gui"></a>
A graphical user interface,
usually controlled by using a mouse.
See also: <a href="#cli">command-line interface</a>.

**home directory**: <a name="home-directory"></a>
The default directory associated with an account on a computer system.
By convention, all of a user's files are stored in or below her home directory.

**immutable**: <a name="immutable"></a>
Unchangeable.
The value of immutable data cannot be altered after it has been created.
See also: (mutable)[#mutable].

**FIXME**: <a name="import"></a>

**in-place operator**: <a name="in-place-operator"></a>
An operator such as `+=` that provides a shorthand notation for
the common case in which the variable being assigned to
is also an operand on the right hand side of the assignment.
For example,
the statement `x += 3` means the same thing as `x = x + 3`.

**FIXME**: <a name="index"></a>

**FIXME**: <a name="infective-license"></a>

**inner loop**: <a name="inner-loop"></a>
A loop that is inside another loop.
See also: (outer loop)[#outer-loop].

**integer**: <a name="integer"></a>
A whole number, such as -12343.
See also: [floating-point number](#float).

**invariant**: <a name="invariant"></a>
An expression whose value doesn't change during the execution of a program,
typically used in an [assertion](#assertion).
See also: (precondition)[#precondition], (postcondition)[#postcondition].

**FIXME**: <a name="loop-body"></a>

**FIXME**: <a name="loop-variable"></a>

**member**: <a name="member"></a>
A variable contained within an [object](#object).

**method**: <a name="method"></a>
A function which is tied to a particular (object)[#object].
Each of an object's methods typically implements one of the things it can do,
or one of the questions it can answer.

**mutable**: <a name="mutable"></a>
Changeable.
The value of mutable data can be updated in place.
See also: (immutable)[#immutable].

**FIXME**: <a name="notional-machine"></a>

**FIXME**: <a name="object"></a>

**FIXME**: <a name="orthogonality"></a>

**outer loop**: <a name="outer-loop"></a>
A loop that contains another loop.
See also: (inner loop)[#inner-loop].

**FIXME**: <a name="parameter"></a>

**parent directory**: <a name="parent-directory"></a>
The directory that "contains" the one in question.
Every directory in a file system except the [root directory](#root-directory) has a parent.
A directory's parent is usually referred to using the shorthand notation `..` (pronounced "dot dot").

**pipe**: <a name="pipe"></a>
A connection from the output of one program to the input of another.
When two or more programs are connected in this way, they are called a "pipeline".

**pipe and filter**: <a name="pipe-and-filter"></a>
A model of programming in which [filters](#filter) that process [streams](#stream) of data
are connected end-to-end.
The pipe and filter model is used extensively in the Unix [shell](#shell).

**FIXME**: <a name="postcondition"></a>

**FIXME**: <a name="precondition"></a>

**FIXME**: <a name="prepared-statement"></a>

**primary key**: <a name="primary-key"></a>
One or more (fields)[#field-database] in a (database table)[#table-database]
whose values are guaranteed to be unique for each (record)[#record-database],
i.e.,
whose values uniquely identify the entry.

**process**: <a name="process"></a>
A running instance of a program,
containing code,
variable values,
open files and network connections,
and so on.
Processes are the "actors" that the [operating system](#operating-system) manages;
it typically runs each process for a few milliseconds at a time
to give the impression that they are executing simultaneously.

**FIXME**: <a name="prompt"></a>

**query**: <a name="query"></a>
A database operation that reads values but does not modify anything.
Queries are expressed in a special-purpose language called (SQL)[#sql].

**FIXME**: <a name="range"></a>

**read-eval-print loop** (REPL): <a name="repl"></a>
a [command-line interface](#cli) that reads a command from the user,
executes it,
prints the result,
and waits for another command.

**record** (in a database): <a name="record-database"></a>
A set of related values making up a single entry in a (database table)[#table-database],
typically shown as a row.
See also: (field)[#field-database].

**FIXME**: <a name="redirection"></a>

**referential integrity**: <a name="referential-integrity"></a>
The internal consistency of values in a database.
If an entry in one table contains a (foreign key)[#foreign-key],
but the corresponding (records)[#record-database] don't exist,
referential integrity has been violated.

**FIXME**: <a name="regression"></a>

**regular expressions** (RE): <a name="regular-expression"></a>
A pattern that specifies a set of character strings.
REs are most often used to find sequences of characters in strings.

**relational database**: <a name="relational-database"></a>
A collection of data organized into (tables)[#table-database].

**relative path**: <a name="relative-path"></a>
A [path](#path) that specifies the location of a file or directory
with respect to the [current working directory](#current-working-directory).
Any path that does not begin with a separator character ("/" or "\") is a relative path.
See also: [absolute path](#absolute-path).

**repository**: <a name="repository"></a>
A central storage area where a (version control)[#version-control] system
stores old (revisions)[#revision] of files
and information about who changed what, when.

**FIXME**: <a name="repository-clone"></a>

**FIXME**: <a name="repository-merge"></a>

**FIXME**: <a name="repository-remote"></a>

**resolve**: <a name="resolve"></a>
To eliminate the (conflicts)[#conflict] between two or more incompatible changes to a file or set of files
being managed by a (version control)[#version-control] system.

**FIXME**: <a name="return-statement"></a>

**FIXME**: <a name="revision"></a>

**FIXME**: <a name="rgb"></a>

**root directory**: <a name="root-directory"></a>
The top-most directory in a [filesystem](#filesystem).
Its name is "/" on Unix (including Linux and Mac OS X) and "\" on Microsoft Windows.

**FIXME**: <a name="sentinel-value"></a>

**FIXME**: <a name="sequence"></a>

**FIXME**: <a name="shape"></a>

**shell**: <a name="shell"></a>
A [command-line interface](#cli)
such as Bash (the Bourne-Again Shell)
or the Microsoft Windows DOS shell
that allows a user to interact with the [operating system](#operating-system).

**FIXME**: <a name="shell-quoting"></a>

**FIXME**: <a name="shell-script"></a>

**FIXME**: <a name="silent-failure"></a>

**slice**: <a name="slice"></a>
A regular subsequence of a larger sequence,
such as the first five elements or every second element.

**SQL** (Structured Query Language): <a name="sql"></a>
A special-purpose language for describing operations on (relational databases)[#relational-database].

**FIXME**: <a name="sql-injection"></a>

**stack frame**: <a name="stack-frame"></a>
A data structure that provides storage for a function's local variables.
Each time a function is called,
a new stack frame is created
and put on the top of the (call stack)[#call-stack].
When the function returns,
the stack frame is discarded.

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

**FIXME**: <a name="sub-directory"></a>

**tab completion**: <a name="tab-completion"></a>
FIXME

**FIXME**: <a name="table-database"></a>

**test-driven development** (TDD): <a name="test-driven-development"></a>
The practice of writing unit tests *before* writing the code they test.

**FIXME**: <a name="test-oracle"></a>

**tuple**: <a name="tuple"></a>
An [immutable](#immutable) (sequence)[#sequence] of values.

**FIXME**: <a name="type"></a>

**FIXME**: <a name="variable"></a>

**version control**: <a name="version-control"></a>
A tool for managing changes to a set of files.
Each set of changes creates a new (revision)[#revision] of the files;
the version control system allows users to recover old revisions reliably,
and helps manage conflicting changes made by different users.

**wildcard**: <a name="wildcard"></a>
A character used in pattern matching.
In the Unix shell,
the wildcard "*" matches zero or more characters,
so that `*.txt` matches all files whose names end in `.txt`.
