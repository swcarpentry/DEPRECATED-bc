---
layout: lesson
root: .
title: Glossary
---
**absolute path**: <a name="absolute-path"></a>
A [path](#path) that refers to a particular location in a file system.
Absolute paths are usually written with respect to the file system's
[root directory](#root-directory),
and begin with either "/" (on Unix) or "\" (on Microsoft Windows).
See also: [relative path](#relative-path).

**additive color model**: <a name="additive-color-model"></a>
A way to represent colors as the sum of contributions from primary colors
such as [red, green, and blue](#rgb).

**aggregation function**: <a name="aggregation-function"></a>
A function such as `sum` or `max` that combines many values to produce a single result.

**alias** (a library): <a name="alias-library"></a>
To give a library a nickname while importing it.

**assertion**: <a name="assertion"></a>
An expression which is supposed to be true at a particular point in a program.
Programmers typically put assertions in their code to check for errors;
if the assertion fails (i.e., if the expression evaluates as false),
the program halts and produces an error message.
See also: [invariant](#invariant), [precondition](#precondition), [postcondition](#postcondition).

**assignment**: <a name="assignment"></a>
To give a value a name by associating a variable with it.

**atomic value**: <a name="atomic-value"></a>
A value that cannot be decomposed into smaller pieces.
For example,
the number 12 is usually considered atomic
(unless we are teaching addition to school children,
in which case we might decompose it into tens and ones).

**call stack**: <a name="call-stack"></a>
A data structure inside a running program that keeps track of active function calls.
Each call's variables are stored in a [stack frame](#stack-frame);
a new stack frame is put on top of the stack for each call,
and discarded when the call is finished.

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

**clone** (a repository): <a name="repository-clone"></a>
To make a copy of a [version control repository](#repository).

**comma-separated values** (CSV): <a name="csv"></a>
A common textual representation for tables
in which the values in each row are separated by commas.

**command-line interface** (CLI): <a name="cli"></a>
An interface based on typing commands,
usually at a [REPL](#repl).
See also: [graphical user interface](#gui).

**comment**: <a name="comment"></a>
A remark in a program that is intended to help human readers understand what is going on,
but is ignored by the computer.
Comments in Python, R, and the Unix shell start with a `#` character and run to the end of the line;
comments in SQL start with `--`,
and other languages have other conventions.

**conditional statement**: <a name="conditional-statement"></a>
A statement in a program that might or might not be executed
depending on whether a test is true or false.

**conflict**: <a name="conflict"></a>
A change made by one user of a [version control system](#version-control)
that is incompatible with changes made by other users.
Helping users [resolve](#resolve) conflicts
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

**data type**: <a name="data-type"></a>
A kind of data value,
such as [integer](#integer) or [character string](#string).

**database manager**: <a name="database-manager"></a>
A program that manages a [relational database](#relational-database).

**default parameter value**: <a name="default-parameter-value"></a>
A value to use for a parameter if nothing is specified explicitly.

**defensive programming**: <a name="defensive-programming"></a>
The practice of writing programs that check their own operation to catch errors as early as possible.

**delimiter**: <a name="delimiter"></a>
A character or characters used to separate individual values,
such as the commas between columns in a [CSV](#csv) file.

**disk block**: <a name="disk-block"></a>
The smallest unit of storage that can be allocated on a computer disk.
Disk blocks are typically 512 bytes in size.

**docstring**: <a name="docstring"></a>
Short for "documentation string",
this refers to textual documentation embedded in Python programs.
Unlike comments,
docstrings are preserved in the running program
and can be examined in interactive sessions.

**documentation**: <a name="documentation"></a>
Human-language text written to explain what software does,
how it works,
or how to use it.

**dotted notation**: <a name="dotted-notation"></a>
A two-part notation used in many programming languages
in which `thing.component` refers to the `component` belonging to `thing`.

**empty string**: <a name="empty-string"></a>
A character string containing no characters,
often thought of as the "zero" of text.

**encapsulation**: <a name="encapsulation"></a>
The practice of hiding something's implementation details
so that the rest of a program can worry about *what* it does
rather than *how* it does it.

**field** (of a database): <a name="field-database"></a>
A set of data values of a particular type,
one for each [record](#record-database) in a [table](#table-database).

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
One or more values in a [database table](#table-database)
that identify a [records](#record-database) in another table.

**function call**: <a name="function-call"></a>
A use of a function in another piece of software.

**function body**: <a name="function-body"></a>
The statements that are executed inside a function.

**function composition**: <a name="function-composition"></a>
The immediate application of one function to the result of another,
such as `f(g(x))`.

**graphical user interface** (GUI): <a name="gui"></a>
A graphical user interface,
usually controlled by using a mouse.
See also: [command-line interface](#cli).

**home directory**: <a name="home-directory"></a>
The default directory associated with an account on a computer system.
By convention, all of a user's files are stored in or below her home directory.

**immutable**: <a name="immutable"></a>
Unchangeable.
The value of immutable data cannot be altered after it has been created.
See also: [mutable](#mutable).

**import**: <a name="import"></a>
To load a library into a program.

**in-place operator**: <a name="in-place-operator"></a>
An operator such as `+=` that provides a shorthand notation for
the common case in which the variable being assigned to
is also an operand on the right hand side of the assignment.
For example,
the statement `x += 3` means the same thing as `x = x + 3`.

**index**: <a name="index"></a>
A subscript that specifies the location of a single value in a collection,
such as a single pixel in an image.

**infective license**: <a name="infective-license"></a>
A license such as the [GPL](http://opensource.org/licenses/GPL-3.0)
that compels people who incorporate material into their own work
to place similar sharing requirements on it.

**inner loop**: <a name="inner-loop"></a>
A loop that is inside another loop.
See also: [outer loop](#outer-loop).

**integer**: <a name="integer"></a>
A whole number, such as -12343.
See also: [floating-point number](#float).

**invariant**: <a name="invariant"></a>
An expression whose value doesn't change during the execution of a program,
typically used in an [assertion](#assertion).
See also: [precondition](#precondition), [postcondition](#postcondition).

**loop variable**: <a name="loop-variable"></a>
The variable that keeps track of the progress of the loop.

**member**: <a name="member"></a>
A variable contained within an [object](#object).

**merge** (a repository): <a name="repository-merge"></a>
To reconcile two sets of change to a [repository](#repository).

**method**: <a name="method"></a>
A function which is tied to a particular [object](#object).
Each of an object's methods typically implements one of the things it can do,
or one of the questions it can answer.

**mutable**: <a name="mutable"></a>
Changeable.
The value of mutable data can be updated in place.
See also: [immutable](#immutable).

**notional machine**: <a name="notional-machine"></a>
An abstraction of a computer used to think about what it can and will do.

**object**: <a name="object"></a>
A particular "chunk" of data associated with specific operations called [methods](#method).

**orthogonal**: <a name="orthogonal"></a>
To have meanings or behaviors that are independent of each other.
If a set of concepts or tools are orthogonal,
they can be combined in any way.

**outer loop**: <a name="outer-loop"></a>
A loop that contains another loop.
See also: [inner loop](#inner-loop).

**parameter**: <a name="parameter"></a>
A value passed into a function,
or a variable named in the function's declaration
that is used to hold such a value.

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

**postcondition**: <a name="postcondition"></a>
A condition that a function (or other block of code) guarantees is true
once it has finished running.
Postconditions are often represented using [assertions](#assertion).

**precondition**: <a name="precondition"></a>
A condition that must be true in order for a function (or other block of code) to run correctly.

**prepared statement**: <a name="prepared-statement"></a>
A template for an [SQL](#sql) query in which some values can be filled in.

**primary key**: <a name="primary-key"></a>
One or more [fields](#field-database) in a [database table](#table-database)
whose values are guaranteed to be unique for each [record](#record-database),
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

**prompt**: <a name="prompt"></a>
A character or characters display by a [REPL](#repl) to show that
it is waiting for its next command.

**query**: <a name="query"></a>
A database operation that reads values but does not modify anything.
Queries are expressed in a special-purpose language called [SQL](#sql).

**quoting** (in the shell): <a name="shell-quoting"></a>
Using quotation marks of various kinds to prevent the shell from interpreting special characters.
For example,
to pass the string `*.txt` to a program,
it is usually necessary to write it as `'*.txt'` (with single quotes)
so that the shell will not try to expand the `*` wildcard.

**read-eval-print loop** (REPL): <a name="repl"></a>
a [command-line interface](#cli) that reads a command from the user,
executes it,
prints the result,
and waits for another command.

**record** (in a database): <a name="record-database"></a>
A set of related values making up a single entry in a [database table](#table-database),
typically shown as a row.
See also: [field](#field-database).

**redirect**: <a name="redirect"></a>
To send a command's output to a file rather than to the screen or another command,
or equivalently to read a command's input from a file.

**referential integrity**: <a name="referential-integrity"></a>
The internal consistency of values in a database.
If an entry in one table contains a [foreign key](#foreign-key),
but the corresponding [records](#record-database) don't exist,
referential integrity has been violated.

**regression**: <a name="regression"></a>
To re-introduce a bug that was once fixed.

**regular expressions** (RE): <a name="regular-expression"></a>
A pattern that specifies a set of character strings.
REs are most often used to find sequences of characters in strings.

**relational database**: <a name="relational-database"></a>
A collection of data organized into [tables](#table-database).

**relative path**: <a name="relative-path"></a>
A [path](#path) that specifies the location of a file or directory
with respect to the [current working directory](#current-working-directory).
Any path that does not begin with a separator character ("/" or "\") is a relative path.
See also: [absolute path](#absolute-path).

**remote repository**: <a name="repository-remote"></a>
A version control [repository](#repository) other than the current one
that the current one is somehow connected to or mirroring.

**repository**: <a name="repository"></a>
A storage area where a [version control](#version-control) system
stores old [revisions](#revision) of files
and information about who changed what, when.

**resolve**: <a name="resolve"></a>
To eliminate the [conflicts](#conflict) between two or more incompatible changes to a file or set of files
being managed by a [version control](#version-control) system.

**return statement**: <a name="return-statement"></a>
A statement that causes a function to stop executing and return a value to its caller immediately.

**revision**: <a name="revision"></a>
A recorded state of a [version control](#version-control) [repository](#repository).

**RGB** (red-green-blue): <a name="rgb"></a>
An [additive model](#additive-color-model)
that represents colors as combinations of red, green, and blue.
Each color's value is typically in the range 0..255
(i.e., a one-byte integer).

**root directory**: <a name="root-directory"></a>
The top-most directory in a [filesystem](#filesystem).
Its name is "/" on Unix (including Linux and Mac OS X) and "\" on Microsoft Windows.

**sentinel value**: <a name="sentinel-value"></a>
A value in a collection that has a special meaning,
such as 999 to mean "age unknown".

**sequence**: <a name="sequence"></a>
An ordered collection of values
whose elements can be specified with a single integer index,
such as a vector.

**shape** (of an array): <a name="shape"></a>
An array's dimensions, represented as a vector.
For example, a 5&times;3 array's shape is `(5,3)`.

**shell**: <a name="shell"></a>
A [command-line interface](#cli)
such as Bash (the Bourne-Again Shell)
or the Microsoft Windows DOS shell
that allows a user to interact with the [operating system](#operating-system).

**shell script**: <a name="shell-script"></a>
A set of [shell](#shell) commands stored in a file for re-use.
A shell script is a program executed by the shell;
the name "script" is used for historical reasons.

**silent failure**: <a name="silent-failure"></a>
Failing without producing any warning messages.
Silent failures are hard to detect and debug.

**slice**: <a name="slice"></a>
A regular subsequence of a larger sequence,
such as the first five elements or every second element.

**SQL** (Structured Query Language): <a name="sql"></a>
A special-purpose language for describing operations on [relational databases](#relational-database).

**SQL injection attack**: <a name="sql-injection-attack"></a>
An attack on a program in which the user's input contains malicious SQL statements.
If this text is copied directly into an SQL statement,
it will be executed in the database.

**stack frame**: <a name="stack-frame"></a>
A data structure that provides storage for a function's local variables.
Each time a function is called,
a new stack frame is created
and put on the top of the [call stack](#call-stack).
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

**stride**: <a name="stride"></a>
The offset between successive elements of a [slice](#slice).

**string**: <a name="string"></a>
Short for "character string",
a [sequence](#sequence) of zero or more characters.

**sub-directory**: <a name="sub-directory"></a>
A directory contained within another directory.

**tab completion**: <a name="tab-completion"></a>
A feature provided by many interactive systems in which
pressing the Tab key triggers automatic completion of the current word or command.

**table** (in a database): <a name="table-database"></a>
A set of data in a [relational database](#relational-database)
organized into a set of [records](#record-database),
each having the same named [fields](#field-database).

**test oracle**: <a name="test-oracle"></a>
A program, device, data set, or human being
against which the results of a test can be compared.

**test-driven development** (TDD): <a name="test-driven-development"></a>
The practice of writing unit tests *before* writing the code they test.

**tuple**: <a name="tuple"></a>
An [immutable](#immutable) [sequence](#sequence) of values.

**variable**: <a name="variable"></a>
A name in a program that is associated with a value or a collection of values.

**version control**: <a name="version-control"></a>
A tool for managing changes to a set of files.
Each set of changes creates a new [revision](#revision) of the files;
the version control system allows users to recover old revisions reliably,
and helps manage conflicting changes made by different users.

**wildcard**: <a name="wildcard"></a>
A character used in pattern matching.
In the Unix shell,
the wildcard "*" matches zero or more characters,
so that `*.txt` matches all files whose names end in `.txt`.
