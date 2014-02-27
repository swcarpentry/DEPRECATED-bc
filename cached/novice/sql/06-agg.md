---
layout: lesson
root: ../..
---

## Aggregation


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Define &quot;aggregation&quot; and give examples of its use.</li>
<li>Write queries that compute aggregated values.</li>
<li>Trace the execution of a query that performs aggregation.</li>
<li>Explain how missing data is handled during aggregation.</li>
</ul>
</div>


<div>
<p>We now want to calculate ranges and averages for our data. We know how to select all of the dates from the <code>Visited</code> table:</p>
</div>


<div class="in">
<pre>%load_ext sqlitemagic</pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select dated from Visited;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>1927-02-08</td></tr>
<tr><td>1927-02-10</td></tr>
<tr><td>1939-01-07</td></tr>
<tr><td>1930-01-12</td></tr>
<tr><td>1930-02-26</td></tr>
<tr><td>None</td></tr>
<tr><td>1932-01-14</td></tr>
<tr><td>1932-03-22</td></tr>
</table></pre>
</div>


<div>
<p>but to combine them, wee must use an <a href="../../gloss.html#aggregation-function">aggregation function</a> such as <code>min</code> or <code>max</code>. Each of these functions takes a set of records as input, and produces a single record as output:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select min(dated) from Visited;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>1927-02-08</td></tr>
</table></pre>
</div>


<div>
<p><img src="img/sql-aggregation.svg" alt="SQL Aggregation" /></p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select max(dated) from Visited;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>1939-01-07</td></tr>
</table></pre>
</div>


<div>
<p><code>min</code> and <code>max</code> are just two of the aggregation functions built into SQL. Three others are <code>avg</code>, <code>count</code>, and <code>sum</code>:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select avg(reading) from Survey where quant=&#39;sal&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>7.20333333333</td></tr>
</table></pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select count(reading) from Survey where quant=&#39;sal&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>9</td></tr>
</table></pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select sum(reading) from Survey where quant=&#39;sal&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>64.83</td></tr>
</table></pre>
</div>


<div>
<p>We used <code>count(reading)</code> here, but we could just as easily have counted <code>quant</code> or any other field in the table, or even used <code>count(*)</code>, since the function doesn't care about the values themselves, just how many values there are.</p>
<p>SQL lets us do several aggregations at once. We can, for example, find the range of sensible salinity measurements:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select min(reading), max(reading) from Survey where quant=&#39;sal&#39; and reading&lt;=1.0;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>0.05</td><td>0.21</td></tr>
</table></pre>
</div>


<div>
<p>We can also combine aggregated results with raw results, although the output might surprise you:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select person, count(*) from Survey where quant=&#39;sal&#39; and reading&lt;=1.0;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>lake</td><td>7</td></tr>
</table></pre>
</div>


<div>
<p>Why does Lake's name appear rather than Roerich's or Dyer's? The answer is that when it has to aggregate a field, but isn't told how to, the database manager chooses an actual value from the input set. It might use the first one processed, the last one, or something else entirely.</p>
<p>Another important fact is that when there are no values to aggregate, aggregation's result is &quot;don't know&quot; rather than zero or some other arbitrary value:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select person, max(reading), sum(reading) from Survey where quant=&#39;missing&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>None</td><td>None</td><td>None</td></tr>
</table></pre>
</div>


<div>
<p>One final important feature of aggregation functions is that they are inconsistent with the rest of SQL in a very useful way. If we add two values, and one of them is null, the result is null. By extension, if we use <code>sum</code> to add all the values in a set, and any of those values are null, the result should also be null. It's much more useful, though, for aggregation functions to ignore null values and only combine those that are non-null. This behavior lets us write our queries as:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select min(dated) from Visited;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>1927-02-08</td></tr>
</table></pre>
</div>


<div>
<p>instead of always having to filter explicitly:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select min(dated) from Visited where dated is not null;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>1927-02-08</td></tr>
</table></pre>
</div>


<div>
<p>Aggregating all records at once doesn't always make sense. For example, suppose Gina suspects that there is a systematic bias in her data, and that some scientists' radiation readings are higher than others. We know that this doesn't work:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant=&#39;rad&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>roe</td><td>8</td><td>6.56</td></tr>
</table></pre>
</div>


<div>
<p>because the database manager selects a single arbitrary scientist's name rather than aggregating separately for each scientist. Since there are only five scientists, she could write five queries of the form:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant=&#39;rad&#39;
and   person=&#39;dyer&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>dyer</td><td>2</td><td>8.81</td></tr>
</table></pre>
</div>


<div>
<p>but this would be tedious, and if she ever had a data set with fifty or five hundred scientists, the chances of her getting all of those queries right is small.</p>
<p>What we need to do is tell the database manager to aggregate the hours for each scientist separately using a <code>group by</code> clause:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select   person, count(reading), round(avg(reading), 2)
from     Survey
where    quant=&#39;rad&#39;
group by person;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>dyer</td><td>2</td><td>8.81</td></tr>
<tr><td>lake</td><td>2</td><td>1.82</td></tr>
<tr><td>pb</td><td>3</td><td>6.66</td></tr>
<tr><td>roe</td><td>1</td><td>11.25</td></tr>
</table></pre>
</div>


<div>
<p><code>group by</code> does exactly what its name implies: groups all the records with the same value for the specified field together so that aggregation can process each batch separately. Since all the records in each batch have the same value for <code>person</code>, it no longer matters that the database manager is picking an arbitrary one to display alongside the aggregated <code>reading</code> values.</p>
</div>


<div>
<p>Just as we can sort by multiple criteria at once, we can also group by multiple criteria. To get the average reading by scientist and quantity measured, for example, we just add another field to the <code>group by</code> clause:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
group by person, quant;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>None</td><td>sal</td><td>1</td><td>0.06</td></tr>
<tr><td>None</td><td>temp</td><td>1</td><td>-26.0</td></tr>
<tr><td>dyer</td><td>rad</td><td>2</td><td>8.81</td></tr>
<tr><td>dyer</td><td>sal</td><td>2</td><td>0.11</td></tr>
<tr><td>lake</td><td>rad</td><td>2</td><td>1.82</td></tr>
<tr><td>lake</td><td>sal</td><td>4</td><td>0.11</td></tr>
<tr><td>lake</td><td>temp</td><td>1</td><td>-16.0</td></tr>
<tr><td>pb</td><td>rad</td><td>3</td><td>6.66</td></tr>
<tr><td>pb</td><td>temp</td><td>2</td><td>-20.0</td></tr>
<tr><td>roe</td><td>rad</td><td>1</td><td>11.25</td></tr>
<tr><td>roe</td><td>sal</td><td>2</td><td>32.05</td></tr>
</table></pre>
</div>


<div>
<p>Note that we have added <code>person</code> to the list of fields displayed, since the results wouldn't make much sense otherwise.</p>
<p>Let's go one step further and remove all the entries where we don't know who took the measurement:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
where    person is not null
group by person, quant
order by person, quant;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>dyer</td><td>rad</td><td>2</td><td>8.81</td></tr>
<tr><td>dyer</td><td>sal</td><td>2</td><td>0.11</td></tr>
<tr><td>lake</td><td>rad</td><td>2</td><td>1.82</td></tr>
<tr><td>lake</td><td>sal</td><td>4</td><td>0.11</td></tr>
<tr><td>lake</td><td>temp</td><td>1</td><td>-16.0</td></tr>
<tr><td>pb</td><td>rad</td><td>3</td><td>6.66</td></tr>
<tr><td>pb</td><td>temp</td><td>2</td><td>-20.0</td></tr>
<tr><td>roe</td><td>rad</td><td>1</td><td>11.25</td></tr>
<tr><td>roe</td><td>sal</td><td>2</td><td>32.05</td></tr>
</table></pre>
</div>


<div>
<p>Looking more closely, this query:</p>
<ol style="list-style-type: decimal">
<li><p>selected records from the <code>Survey</code> table where the <code>person</code> field was not null;</p></li>
<li><p>grouped those records into subsets so that the <code>person</code> and <code>quant</code> values in each subset were the same;</p></li>
<li><p>ordered those subsets first by <code>person</code>, and then within each sub-group by <code>quant</code>; and</p></li>
<li><p>counted the number of records in each subset, calculated the average <code>reading</code> in each, and chose a <code>person</code> and <code>quant</code> value from each (it doesn't matter which ones, since they're all equal).</p></li>
</ol>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>How many temperature readings did Frank Pabodie record, and what was their average value?</p></li>
<li><p>The average of a set of values is the sum of the values divided by the number of values. Does this mean that the <code>avg</code> function returns 2.0 or 3.0 when given the values 1.0, <code>null</code>, and 5.0?</p></li>
<li><p>We want to calculate the difference between each individual radiation reading and the average of all the radiation readings. We write the query:</p>
<pre><code>select reading - avg(reading) from Survey where quant=&#39;rad&#39;;</code></pre>
<p>What does this actually produce, and why?</p></li>
<li><p>The function <code>group_concat(field, separator)</code> concatenates all the values in a field using the specified separator character (or ',' if the separator isn't specified). Use this to produce a one-line list of scientists' names, such as:</p>
<pre><code>William Dyer, Frank Pabodie, Anderson Lake, Valentina Roerich, Frank Danforth</code></pre>
<p>Can you find a way to order the list by surname?</p></li>
</ol>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>An aggregation function combines many values to produce a single new value.</li>
<li>Aggregation functions ignore <code>null</code> values.</li>
<li>Aggregation happens after filtering.</li>
</ul>
</div>
