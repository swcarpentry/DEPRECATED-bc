---
layout: lesson
root: ../..
---

## Selecting Data


<div>
<p>In the late 1920s and early 1930s, William Dyer, Frank Pabodie, and Valentina Roerich led expeditions to the <a href="http://en.wikipedia.org/wiki/Pole_of_inaccessibility">Pole of Inaccessibility</a> in the South Pacific, and then onward to Antarctica. Two years ago, their expeditions were found in a storage locker at Miskatonic University. We have scanned and OCR'd the data they contain, and we now want to store that information in a way that will make search and analysis easy.</p>
<p>We basically have three options: text files, a spreadsheet, or a database. Text files are easiest to create, and work well with version control, but then we would then have to build search and analysis tools ourselves. Spreadsheets are good for doing simple analysis, they don't handle large or complex data sets very well. We would therefore like to put this data in a database, and these lessons will show how to do that.</p>
</div>


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Explain the difference between a table, a record, and a field.</li>
<li>Explain the difference between a database and a database manager.</li>
<li>Write a query to select all values for specific fields from a single table.</li>
</ul>
</div>

### A Few Definitions


<div>
<p>A <a href="../../gloss.html#relational-database">relational database</a> is a way to store and manipulate information that is arranged as <a href="../../gloss.html#table-database">tables</a>. Each table has columns (also known as <a href="../../gloss.html#field-database">fields</a>) which describe the data, and rows (also known as <a href="../../gloss.html#record-database">records</a>) which contain the data.</p>
<p>When we are using a spreadsheet, we put formulas into cells to calculate new values based on old ones. When we are using a database, we send commands (usually called <a href="../../gloss.html#query">queries</a>) to a <a href="../../gloss.html#database-manager">database manager</a>: a program that manipulates the database for us. The database manager does whatever lookups and calculations the query specifies, returning the results in a tabular form that we can then use as a starting point for further queries.</p>
<blockquote>
<p>Every database manager—Oracle, IBM DB2, PostgreSQL, MySQL, Microsoft Access, and SQLite—stores data in a different way, so a database created with one cannot be used directly by another. However, every database manager can import and export data in a variety of formats, so it <em>is</em> possible to move information from one to another.</p>
</blockquote>
<p>Queries are written in a language called <a href="../../gloss.html#sql">SQL</a>, which stands for &quot;Structured Query Language&quot;. SQL provides hundreds of different ways to analyze and recombine data; we will only look at a handful, but that handful accounts for most of what scientists do.</p>
<p>The tables below show the database we will use in our examples:</p>
</div>


<div>
<table>
<tr>
<td valign="top">
<p><strong>Person</strong>: people who took readings.</p>
<table>
  <tr> <th>
ident
</th> <th>
personal
</th> <th>
family
</th> </tr>
  <tr> <td>
dyer
</td> <td>
William
</td> <td>
Dyer
</td> </tr>
  <tr> <td>
pb
</td> <td>
Frank
</td> <td>
Pabodie
</td> </tr>
  <tr> <td>
lake
</td> <td>
Anderson
</td> <td>
Lake
</td> </tr>
  <tr> <td>
roe
</td> <td>
Valentina
</td> <td>
Roerich
</td> </tr>
  <tr> <td>
danforth
</td> <td>
Frank
</td> <td>
Danforth
</td> </tr>
</table>

<p><strong>Site</strong>: locations where readings were taken.</p>
<table>
  <tr> <th>
name
</th> <th>
lat
</th> <th>
long
</th> </tr>
  <tr> <td>
DR-1
</td> <td>
-49.85
</td> <td>
-128.57
</td> </tr>
  <tr> <td>
DR-3
</td> <td>
-47.15
</td> <td>
-126.72
</td> </tr>
  <tr> <td>
MSK-4
</td> <td>
-48.87
</td> <td>
-123.4
</td> </tr>
</table>

<p><strong>Visited</strong>: when readings were taken at specific sites.</p>
<table>
  <tr> <th>
ident
</th> <th>
site
</th> <th>
dated
</th> </tr>
  <tr> <td>
619
</td> <td>
DR-1
</td> <td>
1927-02-08
</td> </tr>
  <tr> <td>
622
</td> <td>
DR-1
</td> <td>
1927-02-10
</td> </tr>
  <tr> <td>
734
</td> <td>
DR-3
</td> <td>
1939-01-07
</td> </tr>
  <tr> <td>
735
</td> <td>
DR-3
</td> <td>
1930-01-12
</td> </tr>
  <tr> <td>
751
</td> <td>
DR-3
</td> <td>
1930-02-26
</td> </tr>
  <tr> <td>
752
</td> <td>
DR-3
</td> <td bgcolor="red">
 
</td> </tr>
  <tr> <td>
837
</td> <td>
MSK-4
</td> <td>
1932-01-14
</td> </tr>
  <tr> <td>
844
</td> <td>
DR-1
</td> <td>
1932-03-22
</td> </tr>
</table>
</td>
<td valign="top">
<p><strong>Survey</strong>: the actual readings.</p>
<table>
  <tr> <th>
taken
</th> <th>
person
</th> <th>
quant
</th> <th>
reading
</th> </tr>
  <tr> <td>
619
</td> <td>
dyer
</td> <td>
rad
</td> <td>
9.82
</td> </tr>
  <tr> <td>
619
</td> <td>
dyer
</td> <td>
sal
</td> <td>
0.13
</td> </tr>
  <tr> <td>
622
</td> <td>
dyer
</td> <td>
rad
</td> <td>
7.8
</td> </tr>
  <tr> <td>
622
</td> <td>
dyer
</td> <td>
sal
</td> <td>
0.09
</td> </tr>
  <tr> <td>
734
</td> <td>
pb
</td> <td>
rad
</td> <td>
8.41
</td> </tr>
  <tr> <td>
734
</td> <td>
lake
</td> <td>
sal
</td> <td>
0.05
</td> </tr>
  <tr> <td>
734
</td> <td>
pb
</td> <td>
temp
</td> <td>
-21.5
</td> </tr>
  <tr> <td>
735
</td> <td>
pb
</td> <td>
rad
</td> <td>
7.22
</td> </tr>
  <tr> <td>
735
</td> <td bgcolor="red">
 
</td> <td>
sal
</td> <td>
0.06
</td> </tr>
  <tr> <td>
735
</td> <td bgcolor="red">
 
</td> <td>
temp
</td> <td>
-26.0
</td> </tr>
  <tr> <td>
751
</td> <td>
pb
</td> <td>
rad
</td> <td>
4.35
</td> </tr>
  <tr> <td>
751
</td> <td>
pb
</td> <td>
temp
</td> <td>
-18.5
</td> </tr>
  <tr> <td>
751
</td> <td>
lake
</td> <td>
sal
</td> <td>
0.1
</td> </tr>
  <tr> <td>
752
</td> <td>
lake
</td> <td>
rad
</td> <td>
2.19
</td> </tr>
  <tr> <td>
752
</td> <td>
lake
</td> <td>
sal
</td> <td>
0.09
</td> </tr>
  <tr> <td>
752
</td> <td>
lake
</td> <td>
temp
</td> <td>
-16.0
</td> </tr>
  <tr> <td>
752
</td> <td>
roe
</td> <td>
sal
</td> <td>
41.6
</td> </tr>
  <tr> <td>
837
</td> <td>
lake
</td> <td>
rad
</td> <td>
1.46
</td> </tr>
  <tr> <td>
837
</td> <td>
lake
</td> <td>
sal
</td> <td>
0.21
</td> </tr>
  <tr> <td>
837
</td> <td>
roe
</td> <td>
sal
</td> <td>
22.5
</td> </tr>
  <tr> <td>
844
</td> <td>
roe
</td> <td>
rad
</td> <td>
11.25
</td> </tr>
</table>
</td>
</tr>
</table>

</div>


<div>
<p>Notice that three entries—one in the <code>Visited</code> table, and two in the <code>Survey</code> table—are shown in red because they don't contain any actual data: we'll return to these missing values <a href="#s:null">later</a>. For now, let's write an SQL query that displays scientists' names. We do this using the SQL command <code>select</code>, giving it the names of the columns we want and the table we want them from. Our query and its output look like this:</p>
</div>


<div class="in">
<pre>%load_ext sqlitemagic</pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select family, personal from Person;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>Dyer</td><td>William</td></tr>
<tr><td>Pabodie</td><td>Frank</td></tr>
<tr><td>Lake</td><td>Anderson</td></tr>
<tr><td>Roerich</td><td>Valentina</td></tr>
<tr><td>Danforth</td><td>Frank</td></tr>
</table></pre>
</div>


<div>
<p>The semi-colon at the end of the query tells the database manager that the query is complete and ready to run. We have written our commands and column names in lower case, and the table name in Title Case, but we don't have to: as the example below shows, SQL is <a href="../../gloss.html#case-insensitive">case insensitive</a>.</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
SeLeCt FaMiLy, PeRsOnAl FrOm PeRsOn;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>Dyer</td><td>William</td></tr>
<tr><td>Pabodie</td><td>Frank</td></tr>
<tr><td>Lake</td><td>Anderson</td></tr>
<tr><td>Roerich</td><td>Valentina</td></tr>
<tr><td>Danforth</td><td>Frank</td></tr>
</table></pre>
</div>


<div>
<p>Whatever casing convention you choose, please be consistent: complex queries are hard enough to read without the extra cognitive load of random capitalization.</p>
</div>


<div>
<p>Going back to our query, it's important to understand that the rows and columns in a database table aren't actually stored in any particular order. They will always be <em>displayed</em> in some order, but we can control that in various ways. For example, we could swap the columns in the output by writing our query as:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select personal, family from Person;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>William</td><td>Dyer</td></tr>
<tr><td>Frank</td><td>Pabodie</td></tr>
<tr><td>Anderson</td><td>Lake</td></tr>
<tr><td>Valentina</td><td>Roerich</td></tr>
<tr><td>Frank</td><td>Danforth</td></tr>
</table></pre>
</div>


<div>
<p>or even repeat columns:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select ident, ident, ident from Person;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>dyer</td><td>dyer</td><td>dyer</td></tr>
<tr><td>pb</td><td>pb</td><td>pb</td></tr>
<tr><td>lake</td><td>lake</td><td>lake</td></tr>
<tr><td>roe</td><td>roe</td><td>roe</td></tr>
<tr><td>danforth</td><td>danforth</td><td>danforth</td></tr>
</table></pre>
</div>


<div>
<p>As a shortcut, we can select all of the columns in a table using <code>*</code>:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Person;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
<tr><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
</table></pre>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Write a query that selects only site names from the <code>Site</code> table.</p></li>
<li><p>Many people format queries as:</p>
<pre><code>SELECT personal, family FROM person;</code></pre>
<p>or as:</p>
<pre><code>select Personal, Family from PERSON;</code></pre>
<p>What style do you find easiest to read, and why?</p></li>
</ol>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>A relational database stores information in tables, each of which has a fixed set of columns and a variable number of records.</li>
<li>A database manager is a program that manipulates information stored in a database.</li>
<li>We write queries in a specialized language called SQL to extract information from databases.</li>
<li>SQL is case-insensitive.</li>
</ul>
</div>
