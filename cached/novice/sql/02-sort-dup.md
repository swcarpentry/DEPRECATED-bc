---
layout: lesson
root: ../..
---

## Sorting and Removing Duplicates


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Write queries that display results in a particular order.</li>
<li>Write queries that eliminate duplicate values from data.</li>
</ul>
</div>


<div>
<p>Data is often redundant, so queries often return redundant information. For example, if we select the quantitites that have been measured from the <code>survey</code> table, we get this:</p>
</div>


<div class="in">
<pre>%load_ext sqlitemagic</pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select quant from Survey;</pre>
</div>

<div class="out">
<pre><table>
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
</table></pre>
</div>


<div>
<p>We can eliminate the redundant output to make the result more readable by adding the <code>distinct</code> keyword to our query:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select distinct quant from Survey;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>rad</td></tr>
<tr><td>sal</td></tr>
<tr><td>temp</td></tr>
</table></pre>
</div>


<div>
<p>If we select more than one column—for example, both the survey site ID and the quantity measured—then the distinct pairs of values are returned:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select distinct taken, quant from Survey;</pre>
</div>

<div class="out">
<pre><table>
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
</table></pre>
</div>


<div>
<p>Notice in both cases that duplicates are removed even if they didn't appear to be adjacent in the database. Again, it's important to remember that rows aren't actually ordered: they're just displayed that way.</p>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li>Write a query that selects distinct dates from the <code>Site</code> table.</li>
</ol>
</div>


<div>
<p>As we mentioned earlier, database records are not stored in any particular order. This means that query results aren't necessarily sorted, and even if they are, we often want to sort them in a different way, e.g., by the name of the project instead of by the name of the scientist. We can do this in SQL by adding an <code>order by</code> clause to our query:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Person order by ident;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
</table></pre>
</div>


<div>
<p>By default, results are sorted in ascending order (i.e., from least to greatest). We can sort in the opposite order using <code>desc</code> (for &quot;descending&quot;):</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from person order by ident desc;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
<tr><td>pb</td><td>Frank</td><td>Pabodie</td></tr>
<tr><td>lake</td><td>Anderson</td><td>Lake</td></tr>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>danforth</td><td>Frank</td><td>Danforth</td></tr>
</table></pre>
</div>


<div>
<p>(And if we want to make it clear that we're sorting in ascending order, we can use <code>asc</code> instead of <code>desc</code>.)</p>
<p>We can also sort on several fields at once. For example, this query sorts results first in ascending order by <code>taken</code>, and then in descending order by <code>person</code> within each group of equal <code>taken</code> values:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select taken, person from Survey order by taken asc, person desc;</pre>
</div>

<div class="out">
<pre><table>
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
</table></pre>
</div>


<div>
<p>This is easier to understand if we also remove duplicates:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select distinct taken, person from Survey order by taken asc, person desc;</pre>
</div>

<div class="out">
<pre><table>
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
</table></pre>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Write a query that returns the distinct dates in the <code>Visited</code> table.</p></li>
<li><p>Write a query that displays the full names of the scientists in the <code>Person</code> table, ordered by family name.</p></li>
</ol>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>The records in a database table are not intrinsically ordered: if we want to display them in some order, we must specify that explicitly.</li>
<li>The values in a database are not guaranteed to be unique: if we want to eliminate duplicates, we must specify that explicitly as well.</li>
</ul>
</div>
