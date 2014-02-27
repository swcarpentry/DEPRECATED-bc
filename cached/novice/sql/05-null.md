---
layout: lesson
root: ../..
---

## Missing Data


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Explain how databases represent missing information.</li>
<li>Explain the three-valued logic databases use when manipulating missing information.</li>
<li>Write queries that handle missing information correctly.</li>
</ul>
</div>


<div>
<p>Real-world data is never complete—there are always holes. Databases represent these holes using special value called <code>null</code>. <code>null</code> is not zero, <code>False</code>, or the empty string; it is a one-of-a-kind value that means &quot;nothing here&quot;. Dealing with <code>null</code> requires a few special tricks and some careful thinking.</p>
<p>To start, let's have a look at the <code>Visited</code> table. There are eight records, but #752 doesn't have a date—or rather, its date is null:</p>
</div>


<div class="in">
<pre>%load_ext sqlitemagic</pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>752</td><td>DR-3</td><td>None</td></tr>
<tr><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></pre>
</div>


<div>
<p>Null doesn't behave like other values. If we select the records that come before 1930:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited where dated&lt;&#39;1930-00-00&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
</table></pre>
</div>


<div>
<p>we get two results, and if we select the ones that come during or after 1930:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited where dated&gt;=&#39;1930-00-00&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></pre>
</div>


<div>
<p>we get five, but record #752 isn't in either set of results. The reason is that <code>null&lt;'1930-00-00'</code> is neither true nor false: null means, &quot;We don't know,&quot; and if we don't know the value on the left side of a comparison, we don't know whether the comparison is true or false. Since databases represent &quot;don't know&quot; as null, the value of <code>null&lt;'1930-00-00'</code> is actually <code>null</code>. <code>null&gt;='1930-00-00'</code> is also null because we can't answer to that question either. And since the only records kept by a <code>where</code> are those for which the test is true, record #752 isn't included in either set of results.</p>
<p>Comparisons aren't the only operations that behave this way with nulls. <code>1+null</code> is <code>null</code>, <code>5*null</code> is <code>null</code>, <code>log(null)</code> is <code>null</code>, and so on. In particular, comparing things to null with = and != produces null:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited where dated=NULL;</pre>
</div>

<div class="out">
<pre><table>

</table></pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited where dated!=NULL;</pre>
</div>

<div class="out">
<pre><table>

</table></pre>
</div>


<div>
<p>To check whether a value is <code>null</code> or not, we must use a special test <code>is null</code>:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited where dated is NULL;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>752</td><td>DR-3</td><td>None</td></tr>
</table></pre>
</div>


<div>
<p>or its inverse <code>is not null</code>:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Visited where dated is not NULL;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></pre>
</div>


<div>
<p>Null values cause headaches wherever they appear. For example, suppose we want to find all the salinity measurements that weren't taken by Dyer. It's natural to write the query like this:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Survey where quant=&#39;sal&#39; and person!=&#39;lake&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>619</td><td>dyer</td><td>sal</td><td>0.13</td></tr>
<tr><td>622</td><td>dyer</td><td>sal</td><td>0.09</td></tr>
<tr><td>752</td><td>roe</td><td>sal</td><td>41.6</td></tr>
<tr><td>837</td><td>roe</td><td>sal</td><td>22.5</td></tr>
</table></pre>
</div>


<div>
<p>but this query filters omits the records where we don't know who took the measurement. Once again, the reason is that when <code>person</code> is <code>null</code>, the <code>!=</code> comparison produces <code>null</code>, so the record isn't kept in our results. If we want to keep these records we need to add an explicit check:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Survey where quant=&#39;sal&#39; and (person!=&#39;lake&#39; or person is null);</pre>
</div>

<div class="out">
<pre><table>
<tr><td>619</td><td>dyer</td><td>sal</td><td>0.13</td></tr>
<tr><td>622</td><td>dyer</td><td>sal</td><td>0.09</td></tr>
<tr><td>735</td><td>None</td><td>sal</td><td>0.06</td></tr>
<tr><td>752</td><td>roe</td><td>sal</td><td>41.6</td></tr>
<tr><td>837</td><td>roe</td><td>sal</td><td>22.5</td></tr>
</table></pre>
</div>


<div>
<p>We still have to decide whether this is the right thing to do or not. If we want to be absolutely sure that we aren't including any measurements by Lake in our results, we need to exclude all the records for which we don't know who did the work.</p>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Write a query that sorts the records in <code>Visited</code> by date, omitting entries for which the date is not known (i.e., is null).</p></li>
<li><p>What do you expect the query:</p>
<pre><code>select * from Visited where dated in (&#39;1927-02-08&#39;, null);</code></pre>
<p>to produce? What does it actually produce?</p></li>
<li><p>Some database designers prefer to use a <a href="../../gloss.html#sentinel-value">sentinel value</a> to mark missing data rather than <code>null</code>. For example, they will use the date &quot;0000-00-00&quot; to mark a missing date, or -1.0 to mark a missing salinity or radiation reading (since actual readings cannot be negative). What does this simplify? What burdens or risks does it introduce?</p></li>
</ol>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>Databases use <code>null</code> to represent missing information.</li>
<li>Any arithmetic or Boolean operation involving <code>null</code> produces <code>null</code> as a result.</li>
<li>The only operators that can safely be used with <code>null</code> are <code>is null</code> and <code>is not null</code>.</li>
</ul>
</div>
