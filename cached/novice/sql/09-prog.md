---
layout: lesson
root: ../..
---

## Programming with Databases


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Write short programs that execute SQL queries.</li>
<li>Trace the execution of a program that contains an SQL query.</li>
<li>Explain why most database applications are written in a general-purpose language rather than in SQL.</li>
</ul>
</div>


<div>
<p>To close, let's have a look at how to access a database from a general-purpose programming language like Python. Other languages use almost exactly the same model: library and function names may differ, but the concepts are the same.</p>
<p>Here's a short Python program that selects latitudes and longitudes from an SQLite database stored in a file called <code>survey.db</code>:</p>
</div>


<div class="in">
<pre>import sqlite3
connection = sqlite3.connect(&#34;survey.db&#34;)
cursor = connection.cursor()
cursor.execute(&#34;select site.lat, site.long from site;&#34;)
results = cursor.fetchall()
for r in results:
    print r
cursor.close()
connection.close()</pre>
</div>

<div class="out">
<pre>(-49.85, -128.57)
(-47.15, -126.72)
(-48.87, -123.4)
</pre>
</div>


<div>
<p>The program starts by importing the <code>sqlite3</code> library. If we were connecting to MySQL, DB2, or some other database, we would import a different library, but all of them provide the same functions, so that the rest of our program does not have to change (at least, not much) if we switch from one database to another.</p>
<p>Line 2 establishes a connection to the database. Since we're using SQLite, all we need to specify is the name of the database file. Other systems may require us to provide a username and password as well. Line 3 then uses this connection to create a <a href="../../gloss.html#cursor">cursor</a>; just like the cursor in an editor, its role is to keep track of where we are in the database.</p>
<p>On line 4, we use that cursor to ask the database to execute a query for us. The query is written in SQL, and passed to <code>cursor.execute</code> as a string. It's our job to make sure that SQL is properly formatted; if it isn't, or if something goes wrong when it is being executed, the database will report an error.</p>
<p>The database returns the results of the query to us in response to the <code>cursor.fetchall</code> call on line 5. This result is a list with one entry for each record in the result set; if we loop over that list (line 6) and print those list entries (line 7), we can see that each one is a tuple with one element for each field we asked for.</p>
<p>Finally, lines 8 and 9 close our cursor and our connection, since the database can only keep a limited number of these open at one time. Since establishing a connection takes time, though, we shouldn't open a connection, do one operation, then close the connection, only to reopen it a few microseconds later to do another operation. Instead, it's normal to create one connection that stays open for the lifetime of the program.</p>
</div>


<div>
<p>Queries in real applications will often depend on values provided by users. For example, this function takes a user's ID as a parameter and returns their name:</p>
</div>


<div class="in">
<pre>def get_name(database_file, person_ident):
    query = &#34;select personal || &#39; &#39; || family from Person where ident=&#39;&#34; + person_ident + &#34;&#39;;&#34;

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query)
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

print &#34;full name for dyer:&#34;, get_name(&#39;survey.db&#39;, &#39;dyer&#39;)</pre>
</div>

<div class="out">
<pre>full name for dyer: William Dyer
</pre>
</div>


<div>
<p>We use string concatenation on the first line of this function to construct a query containing the user ID we have been given. This seems simple enough, but what happens if someone gives us this string as input?</p>
<pre><code>dyer&#39;; drop table Survey; select &#39;</code></pre>
<p>It looks like there's garbage after the name of the project, but it is very carefully chosen garbage. If we insert this string into our query, the result is:</p>
<pre><code>select personal || &#39; &#39; || family from Person where ident=&#39;dyer&#39;; drop tale Survey; select &#39;&#39;;</code></pre>
<p>If we execute this, it will erase one of the tables in our database.</p>
<p>This is called an <a href="../../gloss.html#sql-injection-attack">SQL injection attack</a>, and it has been used to attack thousands of programs over the years. In particular, many web sites that take data from users insert values directly into queries without checking them carefully first.</p>
<p>Since a villain might try to smuggle commands into our queries in many different ways, the safest way to deal with this threat is to replace characters like quotes with their escaped equivalents, so that we can safely put whatever the user gives us inside a string. We can do this by using a <a href="../../gloss.html#prepared-statement">prepared statement</a> instead of formatting our statements as strings. Here's what our example program looks like if we do this:</p>
</div>


<div class="in">
<pre>def get_name(database_file, person_ident):
    query = &#34;select personal || &#39; &#39; || family from Person where ident=?;&#34;

    connection = sqlite3.connect(database_file)
    cursor = connection.cursor()
    cursor.execute(query, [person_ident])
    results = cursor.fetchall()
    cursor.close()
    connection.close()

    return results[0][0]

print &#34;full name for dyer:&#34;, get_name(&#39;survey.db&#39;, &#39;dyer&#39;)</pre>
</div>

<div class="out">
<pre>full name for dyer: William Dyer
</pre>
</div>


<div>
<p>The key changes are in the query string and the <code>execute</code> call. Instead of formatting the query ourselves, we put question marks in the query template where we want to insert values. When we call <code>execute</code>, we provide a list that contains as many values as there are question marks in the query. The library matches values to question marks in order, and translates any special characters in the values into their escaped equivalents so that they are safe to use.</p>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Write a Python program that creates a new database in a file called <code>original.db</code> containing a single table called <code>Pressure</code>, with a single field called <code>reading</code>, and inserts 100,000 random numbers between 10.0 and 25.0. How long does it take this program to run? How long does it take to run a program that simply writes those random numbers to a file?</p></li>
<li><p>Write a Python program that creates a new database called <code>backup.db</code> with the same structure as <code>original.db</code> and copies all the values greater than 20.0 from <code>original.db</code> to <code>backup.db</code>. Which is faster: filtering values in the query, or reading everything into memory and filtering in Python?</p></li>
</ol>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>We usually write database applications in a general-purpose language, and embed SQL queries in it.</li>
<li>To connect to a database, a program must use a library specific to that database manager.</li>
<li>A program may open one or more connections to a single database, and have one or more cursors active in each.</li>
<li>Programs can read query results in batches or all at once.</li>
</ul>
</div>
