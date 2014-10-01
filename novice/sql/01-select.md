---
layout: lesson
root: ../..
---

## Selecting Data


In the late 1920s and early 1930s,
William Dyer,
Frank Pabodie,
and Valentina Roerich led expeditions to the
[Pole of Inaccessibility](http://en.wikipedia.org/wiki/Pole_of_inaccessibility)
in the South Pacific,
and then onward to Antarctica.
Two years ago,
their expeditions were found in a storage locker at Miskatonic University.
We have scanned and OCR'd the data they contain,
and we now want to store that information
in a way that will make search and analysis easy.

We basically have three options:
text files,
a spreadsheet,
or a database.
Text files are easiest to create,
and work well with version control,
but then we would then have to build search and analysis tools ourselves.
Spreadsheets are good for doing simple analysis,
they don't handle large or complex data sets very well.
We would therefore like to put this data in a database,
and these lessons will show how to do that.


<div class="objectives" markdown="1">
#### Objectives

*   Explain the difference between a table, a record, and a field.
*   Explain the difference between a database and a database manager.
*   Write a query to select all values for specific fields from a single table.
</div>

### A Few Definitions


A [relational database](../../gloss.html#relational-database)
is a way to store and manipulate information
that is arranged as [tables](../../gloss.html#table).
Each table has columns (also known as [fields](../../gloss.html#field)) which describe the data,
and rows (also known as [records](../../gloss.html#record)) which contain the data.
  
When we are using a spreadsheet,
we put formulas into cells to calculate new values based on old ones.
When we are using a database,
we send commands
(usually called [queries](../../gloss.html#query))
to a [database manager](../../gloss.html#database-manager):
a program that manipulates the database for us.
The database manager does whatever lookups and calculations the query specifies,
returning the results in a tabular form
that we can then use as a starting point for further queries.
  
> Every database manager&mdash;Oracle,
> IBM DB2, PostgreSQL, MySQL, Microsoft Access, and SQLite&mdash;stores
> data in a different way,
> so a database created with one cannot be used directly by another.
> However,
> every database manager can import and export data in a variety of formats,
> so it *is* possible to move information from one to another.

Queries are written in a language called [SQL](../../gloss.html#sql),
which stands for "Structured Query Language".
SQL provides hundreds of different ways to analyze and recombine data;
we will only look at a handful,
but that handful accounts for most of what scientists do.

The tables below show the database we will use in our examples:


<table>
<tr>
<td valign="top">
<strong>Person</strong>: people who took readings.

<table>
  <tr> <th>ident</th> <th>personal</th> <th>family</th> </tr>
  <tr> <td>dyer</td> <td>William</td> <td>Dyer</td> </tr>
  <tr> <td>pb</td> <td>Frank</td> <td>Pabodie</td> </tr>
  <tr> <td>lake</td> <td>Anderson</td> <td>Lake</td> </tr>
  <tr> <td>roe</td> <td>Valentina</td> <td>Roerich</td> </tr>
  <tr> <td>danforth</td> <td>Frank</td> <td>Danforth</td> </tr>
</table>

<strong>Site</strong>: locations where readings were taken.

<table>
  <tr> <th>name</th> <th>lat</th> <th>long</th> </tr>
  <tr> <td>DR-1</td> <td>-49.85</td> <td>-128.57</td> </tr>
  <tr> <td>DR-3</td> <td>-47.15</td> <td>-126.72</td> </tr>
  <tr> <td>MSK-4</td> <td>-48.87</td> <td>-123.4</td> </tr>
</table>

<strong>Visited</strong>: when readings were taken at specific sites.

<table>
  <tr> <th>ident</th> <th>site</th> <th>dated</th> </tr>
  <tr> <td>619</td> <td>DR-1</td> <td>1927-02-08</td> </tr>
  <tr> <td>622</td> <td>DR-1</td> <td>1927-02-10</td> </tr>
  <tr> <td>734</td> <td>DR-3</td> <td>1939-01-07</td> </tr>
  <tr> <td>735</td> <td>DR-3</td> <td>1930-01-12</td> </tr>
  <tr> <td>751</td> <td>DR-3</td> <td>1930-02-26</td> </tr>
  <tr> <td>752</td> <td>DR-3</td> <td bgcolor="red">&nbsp;</td> </tr>
  <tr> <td>837</td> <td>MSK-4</td> <td>1932-01-14</td> </tr>
  <tr> <td>844</td> <td>DR-1</td> <td>1932-03-22</td> </tr>
</table>
</td>
<td valign="top">
<strong>Survey</strong>: the actual readings.

<table>
  <tr> <th>taken</th> <th>person</th> <th>quant</th> <th>reading</th> </tr>
  <tr> <td>619</td> <td>dyer</td> <td>rad</td> <td>9.82</td> </tr>
  <tr> <td>619</td> <td>dyer</td> <td>sal</td> <td>0.13</td> </tr>
  <tr> <td>622</td> <td>dyer</td> <td>rad</td> <td>7.8</td> </tr>
  <tr> <td>622</td> <td>dyer</td> <td>sal</td> <td>0.09</td> </tr>
  <tr> <td>734</td> <td>pb</td> <td>rad</td> <td>8.41</td> </tr>
  <tr> <td>734</td> <td>lake</td> <td>sal</td> <td>0.05</td> </tr>
  <tr> <td>734</td> <td>pb</td> <td>temp</td> <td>-21.5</td> </tr>
  <tr> <td>735</td> <td>pb</td> <td>rad</td> <td>7.22</td> </tr>
  <tr> <td>735</td> <td bgcolor="red">&nbsp;</td> <td>sal</td> <td>0.06</td> </tr>
  <tr> <td>735</td> <td bgcolor="red">&nbsp;</td> <td>temp</td> <td>-26.0</td> </tr>
  <tr> <td>751</td> <td>pb</td> <td>rad</td> <td>4.35</td> </tr>
  <tr> <td>751</td> <td>pb</td> <td>temp</td> <td>-18.5</td> </tr>
  <tr> <td>751</td> <td>lake</td> <td>sal</td> <td>0.1</td> </tr>
  <tr> <td>752</td> <td>lake</td> <td>rad</td> <td>2.19</td> </tr>
  <tr> <td>752</td> <td>lake</td> <td>sal</td> <td>0.09</td> </tr>
  <tr> <td>752</td> <td>lake</td> <td>temp</td> <td>-16.0</td> </tr>
  <tr> <td>752</td> <td>roe</td> <td>sal</td> <td>41.6</td> </tr>
  <tr> <td>837</td> <td>lake</td> <td>rad</td> <td>1.46</td> </tr>
  <tr> <td>837</td> <td>lake</td> <td>sal</td> <td>0.21</td> </tr>
  <tr> <td>837</td> <td>roe</td> <td>sal</td> <td>22.5</td> </tr>
  <tr> <td>844</td> <td>roe</td> <td>rad</td> <td>11.25</td> </tr>
</table>
</td>
</tr>
</table>


Notice that three entries&mdash;one in the `Visited` table,
and two in the `Survey` table&mdash;are shown in red
because they don't contain any actual data:
we'll return to these missing values [later](#s:null).
For now,
let's write an SQL query that displays scientists' names.
We do this using the SQL command `select`,
giving it the names of the columns we want and the table we want them from.
Our query and its output look like this:


<pre class="in"><code>%load_ext sqlitemagic</code></pre>


<pre class="in"><code>%%sqlite survey.db
select family, personal from Person;</code></pre>

<div class="out"><table>
<tr><td>Dyer</td><td>William</td></tr>
<tr><td>Pabodie</td><td>Frank</td></tr>
<tr><td>Lake</td><td>Anderson</td></tr>
<tr><td>Roerich</td><td>Valentina</td></tr>
<tr><td>Danforth</td><td>Frank</td></tr>
</table></div>


The semi-colon at the end of the query
tells the database manager that the query is complete and ready to run.
We have written our commands and column names in lower case,
and the table name in Title Case,
but we don't have to:
as the example below shows,
SQL is [case insensitive](../../gloss.html#case-insensitive).


<pre class="in"><code>%%sqlite survey.db
SeLeCt FaMiLy, PeRsOnAl FrOm PeRsOn;</code></pre>

<div class="out"><table>
<tr><td>Dyer</td><td>William</td></tr>
<tr><td>Pabodie</td><td>Frank</td></tr>
<tr><td>Lake</td><td>Anderson</td></tr>
<tr><td>Roerich</td><td>Valentina</td></tr>
<tr><td>Danforth</td><td>Frank</td></tr>
</table></div>


Whatever casing convention you choose,
please be consistent:
complex queries are hard enough to read without the extra cognitive load of random capitalization.


Going back to our query,
it's important to understand that
the rows and columns in a database table aren't actually stored in any particular order.
They will always be *displayed* in some order,
but we can control that in various ways.
For example,
we could swap the columns in the output by writing our query as:


<pre class="in"><code>%%sqlite survey.db
select personal, family from Person;</code></pre>

<div class="out"><table>
<tr><td>William</td><td>Dyer</td></tr>
<tr><td>Frank</td><td>Pabodie</td></tr>
<tr><td>Anderson</td><td>Lake</td></tr>
<tr><td>Valentina</td><td>Roerich</td></tr>
<tr><td>Frank</td><td>Danforth</td></tr>
</table></div>


or even repeat columns:


<pre class="in"><code>%%sqlite survey.db
select ident, ident, ident from Person;</code></pre>

<div class="out"><table>
<tr><td>dyer</td><td>dyer</td><td>dyer</td></tr>
<tr><td>pb</td><td>pb</td><td>pb</td></tr>
<tr><td>lake</td><td>lake</td><td>lake</td></tr>
<tr><td>roe</td><td>roe</td><td>roe</td></tr>
<tr><td>danforth</td><td>danforth</td><td>danforth</td></tr>
</table></div>


As a shortcut,
we can select all of the columns in a table using `*`:


<pre class="in"><code>%%sqlite survey.db
select * from Person;</code></pre>

<div class="out"><table>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
<tr><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
</table></div>


#### Challenges

1.  Write a query that selects only site names from the `Site` table.

2.  Many people format queries as:

    ~~~
    SELECT personal, family FROM person;
    ~~~

    or as:

    ~~~
    select Personal, Family from PERSON;
    ~~~

    What style do you find easiest to read, and why?


<div class="keypoints" markdown="1">
#### Key Points

*   A relational database stores information in tables,
    each of which has a fixed set of columns and a variable number of records.
*   A database manager is a program that manipulates information stored in a database.
*   We write queries in a specialized language called SQL to extract information from databases.
*   SQL is case-insensitive.
</div>
