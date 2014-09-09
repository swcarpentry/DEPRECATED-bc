---
layout: lesson
root: ../..
---

## Combining Data


<div class="objectives" markdown="1">
#### Objectives

*   Explain the operation of a query that joins two tables.
*   Explain how to restrict the output of a query containing a join to only include meaningful combinations of values.
*   Write queries that join tables on equal keys.
*   Explain what primary and foreign keys are, and why they are useful.
*   Explain what atomic values are, and why database fields should only contain atomic values.
</div>


In order to submit her data to a web site
that aggregates historical meteorological data,
Gina needs to format it as
latitude, longitude, date, quantity, and reading.
However,
her latitudes and longitudes are in the `Site` table,
while the dates of measurements are in the `Visited` table
and the readings themselves are in the `Survey` table.
She needs to combine these tables somehow.

The SQL command to do this is `join`.
To see how it works,
let's start by joining the `Site` and `Visited` tables:


<pre class="in"><code>%load_ext sqlitemagic</code></pre>


<pre class="in"><code>%%sqlite survey.db
select * from Site join Visited;</code></pre>

<div class="out"><table>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>752</td><td>DR-3</td><td>None</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>752</td><td>DR-3</td><td>None</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>752</td><td>DR-3</td><td>None</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></div>


`join` creates
the [cross product](../../gloss.html#cross-product)
of two tables,
i.e.,
it joins each record of one with each record of the other
to give all possible combinations.
Since there are three records in `Site`
and eight in `Visited`,
the join's output has 24 records.
And since each table has three fields,
the output has six fields.
  
What the join *hasn't* done is
figure out if the records being joined have anything to do with each other.
It has no way of knowing whether they do or not until we tell it how.
To do that,
we add a clause specifying that
we're only interested in combinations that have the same site name:


<pre class="in"><code>%%sqlite survey.db
select * from Site join Visited on Site.name=Visited.site;</code></pre>

<div class="out"><table>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>DR-1</td><td>-49.85</td><td>-128.57</td><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>DR-3</td><td>-47.15</td><td>-126.72</td><td>752</td><td>DR-3</td><td>None</td></tr>
<tr><td>MSK-4</td><td>-48.87</td><td>-123.4</td><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
</table></div>


`on` does the same job as `where`:
it only keeps records that pass some test.
(The difference between the two is that `on` filters records
as they're being created,
while `where` waits until the join is done
and then does the filtering.)
Once we add this to our query,
the database manager throws away records
that combined information about two different sites,
leaving us with just the ones we want.
  
Notice that we used `table.field` to specify field names
in the output of the join.
We do this because tables can have fields with the same name,
and we need to be specific which ones we're talking about.
For example,
if we joined the `person` and `visited` tables,
the result would inherit a field called `ident`
from each of the original tables.

We can now use the same dotted notation
to select the three columns we actually want
out of our join:


<pre class="in"><code>%%sqlite survey.db
select Site.lat, Site.long, Visited.dated
from   Site join Visited
on     Site.name=Visited.site;</code></pre>

<div class="out"><table>
<tr><td>-49.85</td><td>-128.57</td><td>1927-02-08</td></tr>
<tr><td>-49.85</td><td>-128.57</td><td>1927-02-10</td></tr>
<tr><td>-49.85</td><td>-128.57</td><td>1932-03-22</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>None</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-01-12</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-02-26</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1939-01-07</td></tr>
<tr><td>-48.87</td><td>-123.4</td><td>1932-01-14</td></tr>
</table></div>


If joining two tables is good,
joining many tables must be better.
In fact,
we can join any number of tables
simply by adding more `join` clauses to our query,
and more `on` tests to filter out combinations of records
that don't make sense:


<pre class="in"><code>%%sqlite survey.db
select Site.lat, Site.long, Visited.dated, Survey.quant, Survey.reading
from   Site join Visited join Survey
on     Site.name=Visited.site
and    Visited.ident=Survey.taken
and    Visited.dated is not null;</code></pre>

<div class="out"><table>
<tr><td>-49.85</td><td>-128.57</td><td>1927-02-08</td><td>rad</td><td>9.82</td></tr>
<tr><td>-49.85</td><td>-128.57</td><td>1927-02-08</td><td>sal</td><td>0.13</td></tr>
<tr><td>-49.85</td><td>-128.57</td><td>1927-02-10</td><td>rad</td><td>7.8</td></tr>
<tr><td>-49.85</td><td>-128.57</td><td>1927-02-10</td><td>sal</td><td>0.09</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1939-01-07</td><td>rad</td><td>8.41</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1939-01-07</td><td>sal</td><td>0.05</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1939-01-07</td><td>temp</td><td>-21.5</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-01-12</td><td>rad</td><td>7.22</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-01-12</td><td>sal</td><td>0.06</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-01-12</td><td>temp</td><td>-26.0</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-02-26</td><td>rad</td><td>4.35</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-02-26</td><td>sal</td><td>0.1</td></tr>
<tr><td>-47.15</td><td>-126.72</td><td>1930-02-26</td><td>temp</td><td>-18.5</td></tr>
<tr><td>-48.87</td><td>-123.4</td><td>1932-01-14</td><td>rad</td><td>1.46</td></tr>
<tr><td>-48.87</td><td>-123.4</td><td>1932-01-14</td><td>sal</td><td>0.21</td></tr>
<tr><td>-48.87</td><td>-123.4</td><td>1932-01-14</td><td>sal</td><td>22.5</td></tr>
<tr><td>-49.85</td><td>-128.57</td><td>1932-03-22</td><td>rad</td><td>11.25</td></tr>
</table></div>


We can tell which records from `Site`, `Visited`, and `Survey`
correspond with each other
because those tables contain
[primary keys](../../gloss.html#primary-key)
and [foreign keys](../../gloss.html#foreign-key).
A primary key is a value,
or combination of values,
that uniquely identifies each record in a table.
A foreign key is a value (or combination of values) from one table
that identifies a unique record in another table.
Another way of saying this is that
a foreign key is the primary key of one table
that appears in some other table.
In our database,
`Person.ident` is the primary key in the `Person` table,
while `Survey.person` is a foreign key
relating the `Survey` table's entries
to entries in `Person`.

Most database designers believe that
every table should have a well-defined primary key.
They also believe that this key should be separate from the data itself,
so that if we ever need to change the data,
we only need to make one change in one place.
One easy way to do this is
to create an arbitrary, unique ID for each record
as we add it to the database.
This is actually very common:
those IDs have names like "student numbers" and "patient numbers",
and they almost always turn out to have originally been
a unique record identifier in some database system or other.
As the query below demonstrates,
SQLite automatically numbers records as they're added to tables,
and we can use those record numbers in queries:


<pre class="in"><code>%%sqlite survey.db
select rowid, * from Person;</code></pre>

<div class="out"><table>
<tr><td>1</td><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>2</td><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>3</td><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>4</td><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
<tr><td>5</td><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
</table></div>

### Data Hygiene


Now that we have seen how joins work,
we can see why the relational model is so useful
and how best to use it.
The first rule is that every value should be [atomic](../../gloss.html#atomic-value),
i.e.,
not contain parts that we might want to work with separately.
We store personal and family names in separate columns instead of putting the entire name in one column
so that we don't have to use substring operations to get the name's components.
More importantly,
we store the two parts of the name separately because splitting on spaces is unreliable:
just think of a name like "Eloise St. Cyr" or "Jan Mikkel Steubart".

The second rule is that every record should have a unique primary key.
This can be a serial number that has no intrinsic meaning,
one of the values in the record (like the `ident` field in the `Person` table),
or even a combination of values:
the triple `(taken, person, quant)` from the `Survey` table uniquely identifies every measurement.

The third rule is that there should be no redundant information.
For example,
we could get rid of the `Site` table and rewrite the `Visited` table like this:

<table>
  <tr> <td>619</td> <td>-49.85</td> <td>-128.57</td> <td>1927-02-08</td> </tr>
  <tr> <td>622</td> <td>-49.85</td> <td>-128.57</td> <td>1927-02-10</td> </tr>
  <tr> <td>734</td> <td>-47.15</td> <td>-126.72</td> <td>1939-01-07</td> </tr>
  <tr> <td>735</td> <td>-47.15</td> <td>-126.72</td> <td>1930-01-12</td> </tr>
  <tr> <td>751</td> <td>-47.15</td> <td>-126.72</td> <td>1930-02-26</td> </tr>
  <tr> <td>752</td> <td>-47.15</td> <td>-126.72</td> <td>null</td> </tr>
  <tr> <td>837</td> <td>-48.87</td> <td>-123.40</td> <td>1932-01-14</td> </tr>
  <tr> <td>844</td> <td>-49.85</td> <td>-128.57</td> <td>1932-03-22</td> </tr>
</table>

In fact,
we could use a single table that recorded all the information about each reading in each row,
just as a spreadsheet would.
The problem is that it's very hard to keep data organized this way consistent:
if we realize that the date of a particular visit to a particular site is wrong,
we have to change multiple records in the database.
What's worse,
we may have to guess which records to change,
since other sites may also have been visited on that date.

The fourth rule is that the units for every value should be stored explicitly.
Our database doesn't do this,
and that's a problem:
Roerich's salinity measurements are several orders of magnitude larger than anyone else's,
but we don't know if that means she was using parts per million instead of parts per thousand,
or whether there actually was a saline anomaly at that site in 1932.

Stepping back,
data and the tools used to store it have a symbiotic relationship:
we use tables and joins because it's efficient,
provided our data is organized a certain way,
but organize our data that way because we have tools to manipulate it efficiently
if it's in a certain form.
As anthropologists say,
the tool shapes the hand that shapes the tool.


#### Challenges

1.  Write a query that lists all radiation readings from the DR-1 site.

2.  Write a query that lists all sites visited by people named "Frank".

3.  Describe in your own words what the following query produces:

    ~~~
    select Site.name from Site join Visited
    on Site.lat<-49.0 and Site.name=Visited.site and Visited.dated>='1932-00-00';
    ~~~


<div class="keypoints" markdown="1">
#### Key Points

*   Every fact should be represented in a database exactly once.
*   A join produces all combinations of records from one table with records from another.
*   A primary key is a field (or set of fields) whose values uniquely identify the records in a table.
*   A foreign key is a field (or set of fields) in one table whose values are a primary key in another table.
*   We can eliminate meaningless combinations of records by matching primary keys and foreign keys between tables.
*   Keys should be atomic values to make joins simpler and more efficient.
</div>
