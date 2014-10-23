---
layout: lesson
root: .
title: Glossary
---
**<a name="absolute-error">absolute error</a>**:
The absolute value of the difference between a mathematical value
and its finite approximation in a computer.

**<a name="absolute-path">absolute path</a>**:
A [path](#path) that refers to a particular location in a file system.
Absolute paths are usually written with respect to the file system's
[root directory](#root-directory),
and begin with either "/" (on Unix) or "\\" (on Microsoft Windows).
See also: [relative path](#relative-path).

**<a name="access-control-list">access control list</a>** (ACL):
A list of permissions attached to a file or directory
that specifies who can do what with it.

**<a name="additive-color-model">additive color model</a>**:
A way to represent colors as the sum of contributions from primary colors
such as [red, green, and blue](#rgb).

**<a name="aggregation-function">aggregation function</a>**:
A function such as `sum` or `max` that combines many values to produce a single result.

**<a name="alias-library">alias</a>** (a library):
To give a [library](#library) a nickname while importing it.

**<a name="argument">argument</a>**:
A value given to a function or program when it runs.
The term is often used interchangeably (and inconsistently) with [parameter](#parameter).

**<a name="assertion">assertion</a>**:
An expression which is supposed to be true at a particular point in a program.
Programmers typically put assertions in their code to check for errors;
if the assertion fails (i.e., if the expression evaluates as false),
the program halts and produces an error message.
See also: [invariant](#invariant), [precondition](#precondition), [postcondition](#postcondition).

**<a name="assignment">assignment</a>**:
To give a value a name by associating a variable with it.

**<a name="atomic-value">atomic value</a>**:
A value that cannot be decomposed into smaller pieces.
For example,
the number 12 is usually considered atomic
(unless we are teaching addition to school children,
in which case we might decompose it into tens and ones).

**<a name="branch">branch</a>**:
A "parallel universe" in a [version control](#version-control) [repository](#repository).
Programmers typically use branches to isolate different sets of changes from one another during development
so that they can concentrate on one problem at a time.
See also: [merge](#repository-merge).

**<a name="call-stack">call stack</a>**:
A data structure inside a running program that keeps track of active function calls.
Each call's variables are stored in a [stack frame](#stack-frame);
a new stack frame is put on top of the stack for each call,
and discarded when the call is finished.

**<a name="cascading-delete">cascading delete</a>**:
The practice of automatically deleting things in a database
that depend on a record
when that record is deleted.
See also: [referential integrity](#referential-integrity).

**<a name="case-insensitive">case insensitive</a>**:
Treating text as if upper and lower case characters were the same.
See also: [case sensitive](#case-sensitive).

**<a name="catch-exception">catch</a>** (an exception):
To handle an [exception](#exception) that has been [raised](#raise-exception)

**<a name="catch-exception">catch</a>** (an exception):
To handle an [exception](#exception) that has been [raised](#raise-exception)
somewhere else in a program.

**<a name="change-set">change set</a>**:
A group of changes to one or more files
that are [committed](#commit) to a [version control](#version-control) [repository](#repository)
in a single operation.

**<a name="repository-clone">clone</a>** (a repository):
To make a local copy of a [version control repository](#repository).
See also: [fork](#repository-fork).

**<a name="code-review">code review</a>**:
A systematic peer review of a piece of software,
or of changes to a piece of software.
Peer review is often conducted on [pull requests](#pull-request)
before they are [merged](#repository-merge) into a [repository](#repository).

**<a name="csv">comma-separated values</a>** (CSV):
A common textual representation for tables
in which the values in each row are separated by commas.

**<a name="cli">command-line interface</a>** (CLI):
An interface based on typing commands,
usually at a [REPL](#repl).
See also: [graphical user interface](#gui).

**<a name="comment">comment</a>**:
A remark in a program that is intended to help human readers understand what is going on,
but is ignored by the computer.
Comments in Python, R, and the Unix shell start with a `#` character and run to the end of the line;
comments in SQL start with `--`,
and other languages have other conventions.

**<a name="conditional-statement">conditional statement</a>**:
A statement in a program that might or might not be executed
depending on whether a test is true or false.

**<a name="conflict">conflict</a>**:
A change made by one user of a [version control system](#version-control)
that is incompatible with changes made by other users.
Helping users [resolve](#resolve) conflicts
is one of version control's major tasks.

**<a name="cross-product">cross product</a>**:
A pairing of all elements of one set with all elements of another.

**<a name="current-working-directory">current working directory</a>**:
The directory that [relative paths](#relative-path) are calculated from;
equivalently,
the place where files referenced by name only are searched for.
Every [process](#process) has a current working directory.
The current working directory is usually referred to using the shorthand notation `.` (pronounced "dot").

**<a name="cursor">cursor</a>**:
A pointer into a database that keeps track of outstanding operations.

**<a name="data-type">data type</a>**:
A kind of data value,
such as [integer](#integer) or [character string](#string).

**<a name="database-manager">database manager</a>**:
A program that manages a [relational database](#relational-database).

**<a name="default-parameter-value">default parameter value</a>**:
A value to use for a [parameter](#parameter) if nothing is specified explicitly.

**<a name="defensive-programming">defensive programming</a>**:
The practice of writing programs that check their own operation to catch errors as early as possible.

**<a name="delimiter">delimiter</a>**:
A character or characters used to separate individual values,
such as the commas between columns in a [CSV](#csv) file.

**<a name="dependency">dependency</a>**:
Something, usually an input file or software package, that must be available in order to continue a build or install
See [prerequisite](#prerequisite)

**<a name="docstring">docstring</a>**:
Short for "documentation string",
this refers to textual documentation embedded in Python programs.
Unlike comments,
docstrings are preserved in the running program
and can be examined in interactive sessions.

**<a name="documentation">documentation</a>**:
Human-language text written to explain what software does,
how it works,
or how to use it.

**<a name="dotted-notation">dotted notation</a>**:
A two-part notation used in many programming languages
in which `thing.component` refers to the `component` belonging to `thing`.

**<a name="empty-string">empty string</a>**:
A character string containing no characters,
often thought of as the "zero" of text.

**<a name="encapsulation">encapsulation</a>**:
The practice of hiding something's implementation details
so that the rest of a program can worry about *what* it does
rather than *how* it does it.

**<a name="exception">exception</a>**:
An event that disrupts the normal or expected execution of a program.
Most modern languages record information about what went wrong
in a piece of data (also called an exception).
See also: [catch](#catch-exception), [raise](#raise-exception).

**<a name="field-database">field</a>** (of a database):
A set of data values of a particular type,
one for each [record](#record-database) in a [table](#table-database).

**<a name="filename-extension">filename extension</a>**:
The portion of a file's name that comes after the final "." character.
By convention this identifies the file's type:
`.txt` means "text file",
`.png` means "Portable Network Graphics file",
and so on.
These conventions are not enforced by most operating systems:
it is perfectly possible to name an MP3 sound file `homepage.html`.
Since many applications use filename extensions to identify the [MIME type](#mime-type) of the file,
misnaming files may cause those applications to fail.

**<a name="filesystem">filesystem</a>**:
A set of files, directories, and I/O devices (such as keyboards and screens).
A filesystem may be spread across many physical devices,
or many filesystems may be stored on a single physical device;
the [operating system](#operating-system) manages access.

**<a name="filter">filter</a>**:
A program that transforms a stream of data.
Many Unix command-line tools are written as filters:
they read data from [standard input](#standard-input),
process it,
and write the result to [standard output](#standard-output).

**<a name="command-line-flag">flag</a>**:
A terse way to specify an option or setting to a command-line program.
By convention Unix applications use a dash followed by a single letter,
such as `-v`,
or two dashes followed by a word,
such as `--verbose`,
while DOS applications use a slash,
such as `/V`.
Depending on the application, a flag may be followed by a single argument, as in `-o /tmp/output.txt`.

**<a name="float">floating point number</a>** (float):
A number containing a fractional part and an exponent.
See also: [integer](#integer).

**<a name="for-loop">for loop</a>**:
A loop that is executed once for each value in some kind of set, list, or range.
See also: [while loop](#while-loop).

**<a name="foreign-key">foreign key</a>**:
One or more values in a [database table](#table-database)
that identify a [records](#record-database) in another table.

**<a name="repository-fork">fork</a>**:
To [clone](#repository-clone) a [version control](#version-control) [repository](#repository)
on a server.

**<a name="function-body">function body</a>**:
The statements that are executed inside a function.

**<a name="function-call">function call</a>**:
A use of a function in another piece of software.

**<a name="function-composition">function composition</a>**:
The immediate application of one function to the result of another,
such as `f(g(x))`.

**<a name="graph">graph</a>**:
A data structure in software consisting of nodes (containing useful data)
connected by edges (implying relationships between nodes).

**<a name="gui">graphical user interface</a>** (GUI):
A graphical user interface,
usually controlled by using a mouse.
See also: [command-line interface](#cli).

**<a name="home-directory">home directory</a>**:
The default directory associated with an account on a computer system.
By convention, all of a user's files are stored in or below her home directory.

**<a name="http">HTTP</a>**:
The Hypertext Transfer [Protocol](#protocol) used for sharing web pages and other data
on the World Wide Web.

**<a name="immutable">immutable</a>**:
Unchangeable.
The value of immutable data cannot be altered after it has been created.
See also: [mutable](#mutable).

**<a name="import">import</a>**:
To load a [library](#library) into a program.

**<a name="in-place-operator">in-place operator</a>**:
An operator such as `+=` that provides a shorthand notation for
the common case in which the variable being assigned to
is also an operand on the right hand side of the assignment.
For example,
the statement `x += 3` means the same thing as `x = x + 3`.

**<a name="index">index</a>**:
A subscript that specifies the location of a single value in a collection,
such as a single pixel in an image.

**<a name="infective-license">infective license</a>**:
A license such as the [GPL](http://opensource.org/licenses/GPL-3.0)
that compels people who incorporate material into their own work
to place similar sharing requirements on it.

**<a name="inner-loop">inner loop</a>**:
A loop that is inside another loop.
See also: [outer loop](#outer-loop).

**<a name="integer">integer</a>**:
A whole number, such as -12343.
See also: [floating-point number](#float).

**<a name="invariant">invariant</a>**:
An expression whose value doesn't change during the execution of a program,
typically used in an [assertion](#assertion).
See also: [precondition](#precondition), [postcondition](#postcondition).

**<a name="library">library</a>**:
A family of code units (functions, classes, variables) that implement a set of
related tasks.

**<a name="loop-body">loop body</a>**:
The set of statements or commands that are repeated inside a [for loop](#for-loop)
or [while loop](#while-loop).

**<a name="loop-variable">loop variable</a>**:
The variable that keeps track of the progress of the loop.

**<a name="macro">macro</a>**:
In [Makefiles](#makefile), [variables](#variable) are often referred to as macros because they are expanded when read.

**<a name="makefile">makefile</a>**:
An input file to the `make` program. It tells `make` what to do.

**<a name="member">member</a>**:
A variable contained within an [object](#object).

**<a name="repository-merge">merge</a>** (a repository):
To reconcile two sets of change to a [repository](#repository).

**<a name="method">method</a>**:
A function which is tied to a particular [object](#object).
Each of an object's methods typically implements one of the things it can do,
or one of the questions it can answer.

**<a name="mutable">mutable</a>**:
Changeable.
The value of mutable data can be updated in place.
See also: [immutable](#immutable).

**<a name="notional-machine">notional machine</a>**:
An abstraction of a computer used to think about what it can and will do.

**<a name="orthogonal">orthogonal</a>**:
To have meanings or behaviors that are independent of each other.
If a set of concepts or tools are orthogonal,
they can be combined in any way.

**<a name="outer-loop">outer loop</a>**:
A loop that contains another loop.
See also: [inner loop](#inner-loop).

**<a name="parameter">parameter</a>**:
A variable named in the function's declaration that is used to hold a value passed into the call.
The term is often used interchangeably (and inconsistently) with [argument](#argument).

**<a name="parent-directory">parent directory</a>**:
The directory that "contains" the one in question.
Every directory in a file system except the [root directory](#root-directory) has a parent.
A directory's parent is usually referred to using the shorthand notation `..` (pronounced "dot dot").

**<a name="phony-target">phony target</a>**:
A type of [target](#target) within a [Makefile](#makefile) [rule](#make-rule) that does not correspond
to an actual file that needs to be created, but instead serves as a place holder to ensure
[dependencies](#dependencies) get run

**<a name="pipe">pipe</a>**:
A connection from the output of one program to the input of another.
When two or more programs are connected in this way, they are called a "pipeline".

**<a name="pipe-and-filter">pipe and filter</a>**:
A model of programming in which [filters](#filter) that process [streams](#stream) of data
are connected end-to-end.
The pipe and filter model is used extensively in the Unix [shell](#shell).

**<a name="postcondition">postcondition</a>**:
A condition that a function (or other block of code) guarantees is true
once it has finished running.
Postconditions are often represented using [assertions](#assertion).

**<a name="precondition">precondition</a>**:
A condition that must be true in order for a function (or other block of code) to run correctly.

**<a name="prepared-statement">prepared statement</a>**:
A template for an [SQL](#sql) query in which some values can be filled in.

**<a name="prerequisite">prerequisite</a>**:
An input file (or list of input files) defined in a [rule](#make-rule) required to create a [target](#make-target)
See [dependency](#dependency)

**<a name="primary-key">primary key</a>**:
One or more [fields](#field-database) in a [database table](#table-database)
whose values are guaranteed to be unique for each [record](#record-database),
i.e.,
whose values uniquely identify the entry.

**<a name="process">process</a>**:
A running instance of a program,
containing code,
variable values,
open files and network connections,
and so on.
Processes are the "actors" that the [operating system](#operating-system) manages;
it typically runs each process for a few milliseconds at a time
to give the impression that they are executing simultaneously.

**<a name="prompt">prompt</a>**:
A character or characters display by a [REPL](#repl) to show that
it is waiting for its next command.

**<a name="protocol">protocol</a>**:
A set of rules that define how one computer communicates with another.
Common protocols on the Internet include [HTTP](#http) and [SSH](#ssh).

**<a name="pull-request">pull request</a>**:
A set of changes created in one [version control](#version-control) [repository](#repository)
that is being offered to another for [merging](#repository-merge).

**<a name="query">query</a>**:
A database operation that reads values but does not modify anything.
Queries are expressed in a special-purpose language called [SQL](#sql).

**<a name="shell-quoting">quoting</a>** (in the shell):
Using quotation marks of various kinds to prevent the shell from interpreting special characters.
For example,
to pass the string `*.txt` to a program,
it is usually necessary to write it as `'*.txt'` (with single quotes)
so that the shell will not try to expand the `*` wildcard.

**<a name="raise-exception">raise</a>** (an exception):
To explicitly signal that an [exception](#exception) has occured in a program.
See also: [catch](#catch-exception).

**<a name="repl">read-eval-print loop</a>** (REPL):
A [command-line interface](#cli) that reads a command from the user,
executes it,
prints the result,
and waits for another command.

**<a name="recipe">recipe</a>**:
In a [Makefile](#makefile), a recipe is part of a [rule](#make-rule) that describes
programs to run in order to create the [target](#make-target) from the [prerequisites](#prerequisite)

**<a name="record-database">record</a>** (in a database):
A set of related values making up a single entry in a [database table](#table-database),
typically shown as a row.
See also: [field](#field-database).

**<a name="redirect">redirect</a>**:
To send a command's output to a file rather than to the screen or another command,
or equivalently to read a command's input from a file.

**<a name="referential-integrity">referential integrity</a>**:
The internal consistency of values in a database.
If an entry in one table contains a [foreign key](#foreign-key),
but the corresponding [records](#record-database) don't exist,
referential integrity has been violated.

**<a name="regression">regression</a>**:
To re-introduce a bug that was once fixed.

**<a name="regular-expression">regular expressions</a>** (RE):
A pattern that specifies a set of character strings.
REs are most often used to find sequences of characters in strings.

**<a name="relational-database">relational database</a>**:
A collection of data organized into [tables](#table-database).

**<a name="relative-error">relative error</a>**:
The ratio of the [absolute error](#absolute-error) in an approximation of a value
to the value being approximated.

**<a name="relative-path">relative path</a>**:
A [path](#path) that specifies the location of a file or directory
with respect to the [current working directory](#current-working-directory).
Any path that does not begin with a separator character ("/" or "\\") is a relative path.
See also: [absolute path](#absolute-path).

**<a name="remote-login">remote login</a>**:
To connect to a computer over a network,
e.g., to run a [shell](#shell) on it.
See also: [SSH](#ssh).

**<a name="repository-remote">remote repository</a>**:
A version control [repository](#repository) other than the current one
that the current one is somehow connected to or mirroring.

**<a name="repository">repository</a>**:
A storage area where a [version control](#version-control) system
stores old [revisions](#revision) of files
and information about who changed what, when.

**<a name="resolve">resolve</a>**:
To eliminate the [conflicts](#conflict) between two or more incompatible changes to a file or set of files
being managed by a [version control](#version-control) system.

**<a name="return-statement">return statement</a>**:
A statement that causes a function to stop executing and return a value to its caller immediately.

**<a name="revision">revision</a>**:
A recorded state of a [version control](#version-control) [repository](#repository).

**<a name="rgb">RGB</a>** (red-green-blue):
An [additive model](#additive-color-model)
that represents colors as combinations of red, green, and blue.
Each color's value is typically in the range 0..255
(i.e., a one-byte integer).

**<a name="root-directory">root directory</a>**:
The top-most directory in a [filesystem](#filesystem).
Its name is "/" on Unix (including Linux and Mac OS X) and "\\" on Microsoft Windows.

**<a name="make-rule">rule</a>**:
A series of related lines in a [Makefile](#makefile) that define a [recipe](#recipe)
to turn a set of [prerequisites](#prerequisite) into a desired [target](#make-target)

**<a name="search-path">search path</a>**:
The list of directories in which the [shell](#shell) searches for programs when they are run.

**<a name="sentinel-value">sentinel value</a>**:
A value in a collection that has a special meaning,
such as 999 to mean "age unknown".

**<a name="shape">shape</a>** (of an array):
An array's dimensions, represented as a vector.
For example, a 5&times;3 array's shape is `(5,3)`.

**<a name="shell">shell</a>**:
A [command-line interface](#cli)
such as Bash (the Bourne-Again Shell)
or the Microsoft Windows DOS shell
that allows a user to interact with the [operating system](#operating-system).

**<a name="shell-script">shell script</a>**:
A set of [shell](#shell) commands stored in a file for re-use.
A shell script is a program executed by the shell;
the name "script" is used for historical reasons.

**<a name="sign-and-magnitude">sign and magnitude</a>**:
A scheme for representing numbers in which one bit indicates the sign (positive or negative)
and the other bits store the number's absolute value.
See also: [two's complement](#twos-complement).

**<a name="silent-failure">silent failure</a>**:
Failing without producing any warning messages.
Silent failures are hard to detect and debug.

**<a name="slice">slice</a>**:
A regular subsequence of a larger sequence,
such as the first five elements or every second element.

**<a name="sql">SQL</a>** (Structured Query Language):
A special-purpose language for describing operations on [relational databases](#relational-database).

**<a name="sql-injection-attack">SQL injection attack</a>**:
An attack on a program in which the user's input contains malicious SQL statements.
If this text is copied directly into an SQL statement,
it will be executed in the database.

**<a name="ssh">SSH</a>**:
The Secure Shell [protocol](#protocol) used for secure communication between computers.
SSH is often used for [remote login](#remote-login) between computers.

**<a name="ssh-key">SSH key</a>**:
A digital key that identifies one computer or user to another.

**<a name="stack-frame">stack frame</a>**:
A data structure that provides storage for a function's local variables.
Each time a function is called,
a new stack frame is created
and put on the top of the [call stack](#call-stack).
When the function returns,
the stack frame is discarded.

**<a name="standard-input">standard input</a>** (stdin):
A process's default input stream.
In interactive command-line applications,
it is typically connected to the keyboard;;
in a [pipe](#pipe),
it receives data from the [standard output](#standard-output) of the preceding process.

**<a name="standard-output">standard output</a>** (stdout):
A process's default output stream.
In interactive command-line applications,
data sent to standard output is displayed on the screen;
in a [pipe](#pipe),
it is passed to the [standard input](#standard-input) of the next process.

**<a name="stride">stride</a>**:
The offset between successive elements of a [slice](#slice).

**<a name="string">string</a>**:
Short for "character string",
a [sequence](#sequence) of zero or more characters.

**<a name="sub-directory">sub-directory</a>**:
A directory contained within another directory.

**<a name="tab-completion">tab completion</a>**:
A feature provided by many interactive systems in which
pressing the Tab key triggers automatic completion of the current word or command.

**<a name="table-database">table</a>** (in a database):
A set of data in a [relational database](#relational-database)
organized into a set of [records](#record-database),
each having the same named [fields](#field-database).

**<a name="make-target">target</a>**:
The goal of a [rule](#make-rule) in a [Makefile](#makefile).

**<a name="test-oracle">test oracle</a>**:
A program, device, data set, or human being
against which the results of a test can be compared.

**<a name="test-driven-development">test-driven development</a>** (TDD):
The practice of writing unit tests *before* writing the code they test.

**<a name="timestamp">timestamp</a>**:
A record of when a particular event occurred.

**<a name="tuple">tuple</a>**:
An [immutable](#immutable) [sequence](#sequence) of values.

**<a name="twos-complement">two's complement</a>**:
A scheme for representing numbers which wraps around like an odometer
so that 111...111 represents -1.
See also: [sign and magnitude](#sign-and-magnitude).

**<a name="user-group">user group</a>**:
A set of users on a computer system.

**<a name="user-group-id">user group ID</a>**:
A numerical ID that specifies a [user group](#user-group).

**<a name="user-group-name">user group name</a>**:
A textual name for a [user group](#user-group).

**<a name="user-id">user ID</a>**:
A numerical ID that specifies an individual user on a computer system.
See also: [user name](#user-name).

**<a name="user-name">user name</a>**:
A textual name for a user on a computer system.
See also: [user ID](#user-id).

**<a name="variable">variable</a>**:
A name in a program that is associated with a value or a collection of values.

**<a name="version-control">version control</a>**:
A tool for managing changes to a set of files.
Each set of changes creates a new [revision](#revision) of the files;
the version control system allows users to recover old revisions reliably,
and helps manage conflicting changes made by different users.

**<a name="while-loop">while loop</a>**:
A loop that keeps executing as long as some condition is true.
See also: [for loop](#for-loop).

**<a name="wildcard">wildcard</a>**:
A character used in pattern matching.
In the Unix shell,
the wildcard "*" matches zero or more characters,
so that `*.txt` matches all files whose names end in `.txt`.
