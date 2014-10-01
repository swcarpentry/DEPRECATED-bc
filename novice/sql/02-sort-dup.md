---
layout: lesson
root: ../..
---

## Sorting and Removing Duplicates


<div class="objectives" markdown="1">
#### Objectives

*   Write queries that display results in a particular order.
*   Write queries that eliminate duplicate values from data.
</div>


Data is often redundant,
so queries often return redundant information.
For example,
if we select the quantitites that have been measured
from the `survey` table,
we get this:


<pre class="in"><code>%load_ext sqlitemagic</code></pre>


<pre class="in"><code>%%sqlite survey.db
select quant from Survey;</code></pre>

<div class="out"><table>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>temp</td></tr>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>temp</td></tr>
<tr><td>rad</td></tr>
<tr><td>temp</td></tr>
<tr><td>sal</td></tr>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>temp</td></tr>
<tr><td>sal</td></tr>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>sal</td></tr>
<tr><td>rad</td></tr>
</table></div>


We can eliminate the redundant output
to make the result more readable
by adding the `distinct` keyword
to our query:


<pre class="in"><code>%%sqlite survey.db
select distinct quant from Survey;</code></pre>

<div class="out"><table>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>temp</td></tr>
</table></div>


If we select more than one column&mdash;for example,
both the survey site ID and the quantity measured&mdash;then
the distinct pairs of values are returned:


<pre class="in"><code>%%sqlite survey.db
select distinct taken, quant from Survey;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>rad</td></tr>
<tr><td>619</td><td>sal</td></tr>
<tr><td>622</td><td>rad</td></tr>
<tr><td>622</td><td>sal</td></tr>
<tr><td>734</td><td>rad</td></tr>
<tr><td>734</td><td>sal</td></tr>
<tr><td>734</td><td>temp</td></tr>
<tr><td>735</td><td>rad</td></tr>
<tr><td>735</td><td>sal</td></tr>
<tr><td>735</td><td>temp</td></tr>
<tr><td>751</td><td>rad</td></tr>
<tr><td>751</td><td>temp</td></tr>
<tr><td>751</td><td>sal</td></tr>
<tr><td>752</td><td>rad</td></tr>
<tr><td>752</td><td>sal</td></tr>
<tr><td>752</td><td>temp</td></tr>
<tr><td>837</td><td>rad</td></tr>
<tr><td>837</td><td>sal</td></tr>
<tr><td>844</td><td>rad</td></tr>
</table></div>


Notice in both cases that duplicates are removed
even if they didn't appear to be adjacent in the database.
Again,
it's important to remember that rows aren't actually ordered:
they're just displayed that way.


#### Challenges

1.  Write a query that selects distinct dates from the `Site` table.


As we mentioned earlier,
database records are not stored in any particular order.
This means that query results aren't necessarily sorted,
and even if they are,
we often want to sort them in a different way,
e.g., by the name of the project instead of by the name of the scientist.
We can do this in SQL by adding an `order by` clause to our query:


<pre class="in"><code>%%sqlite survey.db
select * from Person order by ident;</code></pre>

<div class="out"><table>
<tr><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
</table></div>


By default,
results are sorted in ascending order
(i.e.,
from least to greatest).
We can sort in the opposite order using `desc` (for "descending"):


<pre class="in"><code>%%sqlite survey.db
select * from person order by ident desc;</code></pre>

<div class="out"><table>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
<tr><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
</table></div>


(And if we want to make it clear that we're sorting in ascending order,
we can use `asc` instead of `desc`.)
  
We can also sort on several fields at once.
For example,
this query sorts results first in ascending order by `taken`,
and then in descending order by `person`
within each group of equal `taken` values:


<pre class="in"><code>%%sqlite survey.db
select taken, person from Survey order by taken asc, person desc;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>dyer</td></tr>
<tr><td>619</td><td>dyer</td></tr>
<tr><td>622</td><td>dyer</td></tr>
<tr><td>622</td><td>dyer</td></tr>
<tr><td>734</td><td>pb</td></tr>
<tr><td>734</td><td>pb</td></tr>
<tr><td>734</td><td>lake</td></tr>
<tr><td>735</td><td>pb</td></tr>
<tr><td>735</td><td>None</td></tr>
<tr><td>735</td><td>None</td></tr>
<tr><td>751</td><td>pb</td></tr>
<tr><td>751</td><td>pb</td></tr>
<tr><td>751</td><td>lake</td></tr>
<tr><td>752</td><td>roe</td></tr>
<tr><td>752</td><td>lake</td></tr>
<tr><td>752</td><td>lake</td></tr>
<tr><td>752</td><td>lake</td></tr>
<tr><td>837</td><td>roe</td></tr>
<tr><td>837</td><td>lake</td></tr>
<tr><td>837</td><td>lake</td></tr>
<tr><td>844</td><td>roe</td></tr>
</table></div>


This is easier to understand if we also remove duplicates:


<pre class="in"><code>%%sqlite survey.db
select distinct taken, person from Survey order by taken asc, person desc;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>dyer</td></tr>
<tr><td>622</td><td>dyer</td></tr>
<tr><td>734</td><td>pb</td></tr>
<tr><td>734</td><td>lake</td></tr>
<tr><td>735</td><td>pb</td></tr>
<tr><td>735</td><td>None</td></tr>
<tr><td>751</td><td>pb</td></tr>
<tr><td>751</td><td>lake</td></tr>
<tr><td>752</td><td>roe</td></tr>
<tr><td>752</td><td>lake</td></tr>
<tr><td>837</td><td>roe</td></tr>
<tr><td>837</td><td>lake</td></tr>
<tr><td>844</td><td>roe</td></tr>
</table></div>


#### Challenges

1.  Write a query that returns the distinct dates in the `Visited` table.

2.  Write a query that displays the full names of the scientists in the `Person` table, ordered by family name.


<div class="keypoints" markdown="1">
#### Key Points

*   The records in a database table are not intrinsically ordered:
    if we want to display them in some order,
    we must specify that explicitly.
*   The values in a database are not guaranteed to be unique:
    if we want to eliminate duplicates,
    we must specify that explicitly as well.
</div>
