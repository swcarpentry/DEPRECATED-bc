---
layout: lesson
root: ../..
---

## Creating and Modifying Data


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Write queries that creates tables.</li>
<li>Write queries to insert, modify, and delete records.</li>
</ul>
</div>


<div>
<p>So far we have only looked at how to get information out of a database, both because that is more frequent than adding information, and because most other operations only make sense once queries are understood. If we want to create and modify data, we need to know two other pairs of commands.</p>
<p>The first pair are <code>create table</code> and <code>drop table</code>. While they are written as two words, they are actually single commands. The first one creates a new table; its arguments are the names and types of the table's columns. For example, the following statements create the four tables in our survey database:</p>
<pre><code>create table Person(ident text, personal text, family text);
create table Site(name text, lat real, long real);
create table Visited(ident integer, site text, dated text);
create table Survey(taken integer, person text, quant real, reading real);</code></pre>
<p>We can get rid of one of our tables using:</p>
<pre><code>drop table Survey;</code></pre>
<p>Be very careful when doing this: most databases have some support for undoing changes, but it's better not to have to rely on it.</p>
<p>Different database systems support different data types for table columns, but most provide the following:</p>
<table>
  <tr> <td>
integer
</td> <td>
a signed integer
</td> </tr>
  <tr> <td>
real
</td> <td>
a floating point number
</td> </tr>
  <tr> <td>
text
</td> <td>
a character string
</td> </tr>
  <tr> <td>
blob
</td> <td>
a &quot;binary large object&quot;, such as an image
</td> </tr>
</table>

<p>Most databases also support Booleans and date/time values; SQLite uses the integers 0 and 1 for the former, and represents the latter as discussed <a href="#a:dates">earlier</a>. An increasing number of databases also support geographic data types, such as latitude and longitude. Keeping track of what particular systems do or do not offer, and what names they give different data types, is an unending portability headache.</p>
<p>When we create a table, we can specify several kinds of constraints on its columns. For example, a better definition for the <code>Survey</code> table would be:</p>
<pre><code>create table Survey(
    taken   integer not null, -- where reading taken
    person  text,             -- may not know who took it
    quant   real not null,    -- the quantity measured
    reading real not null,    -- the actual reading
    primary key(taken, quant),
    foreign key(taken) references Visited(ident),
    foreign key(person) references Person(ident)
);</code></pre>
<p>Once again, exactly what constraints are avialable and what they're called depends on which database manager we are using.</p>
<p>Once tables have been created, we can add and remove records using our other pair of commands, <code>insert</code> and <code>delete</code>. The simplest form of <code>insert</code> statement lists values in order:</p>
<pre><code>insert into Site values(&#39;DR-1&#39;, -49.85, -128.57);
insert into Site values(&#39;DR-3&#39;, -47.15, -126.72);
insert into Site values(&#39;MSK-4&#39;, -48.87, -123.40);</code></pre>
<p>We can also insert values into one table directly from another:</p>
<pre><code>create table JustLatLong(lat text, long text);
insert into JustLatLong select lat, long from site;</code></pre>
<p>Deleting records can be a bit trickier, because we have to ensure that the database remains internally consistent. If all we care about is a single table, we can use the <code>delete</code> command with a <code>where</code> clause that matches the records we want to discard. For example, once we realize that Frank Danforth didn't take any measurements, we can remove him from the <code>Person</code> table like this:</p>
<pre><code>delete from Person where ident = &quot;danforth&quot;;</code></pre>
<p>But what if we removed Anderson Lake instead? Our <code>Survey</code> table would still contain seven records of measurements he'd taken, but that's never supposed to happen: <code>Survey.person</code> is a foreign key into the <code>Person</code> table, and all our queries assume there will be a row in the latter matching every value in the former.</p>
<p>This problem is called <a href="../../gloss.html#referential-integrity">referential integrity</a>: we need to ensure that all references between tables can always be resolved correctly. One way to do this is to delete all the records that use <code>'lake'</code> as a foreign key before deleting the record that uses it as a primary key. If our database manager supports it, we can automate this using <a href="../../gloss.html#cascading-delete">cascading delete</a>. However, this technique is outside the scope of this chapter.</p>
<blockquote>
<p>Many applications use a hybrid storage model instead of putting everything into a database: the actual data (such as astronomical images) is stored in files, while the database stores the files' names, their modification dates, the region of the sky they cover, their spectral characteristics, and so on. This is also how most music player software is built: the database inside the application keeps track of the MP3 files, but the files themselves live on disk.</p>
</blockquote>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Write an SQL statement to replace all uses of <code>null</code> in <code>Survey.person</code> with the string <code>'unknown'</code>.</p></li>
<li><p>One of our colleagues has sent us a <a href="../../gloss.html#csv">CSV</a> file containing temperature readings by Robert Olmstead, which is formatted like this:</p>
<pre><code>Taken,Temp
619,-21.5
622,-15.5</code></pre>
<p>Write a small Python program that reads this file in and prints out the SQL <code>insert</code> statements needed to add these records to the survey database. Note: you will need to add an entry for Olmstead to the <code>Person</code> table. If you are testing your program repeatedly, you may want to investigate SQL's <code>insert or replace</code> command.</p></li>
<li><p>SQLite has several administrative commands that aren't part of the SQL standard. One of them is <code>.dump</code>, which prints the SQL commands needed to re-create the database. Another is <code>.load</code>, which reads a file created by <code>.dump</code> and restores the database. A colleague of yours thinks that storing dump files (which are text) in version control is a good way to track and manage changes to the database. What are the pros and cons of this approach? (Hint: records aren't stored in any particular order.)</p></li>
</ol>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>Database tables are created using queries that specify their names and the names and properties of their fields.</li>
<li>Records can be inserted, updated, or deleted using queries.</li>
<li>It is simpler and safer to modify data when every record has a unique primary key.</li>
</ul>
</div>
