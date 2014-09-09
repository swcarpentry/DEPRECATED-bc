---
layout: lesson
root: ../..
---

## Aggregation


<div class="objectives" markdown="1">
#### Objectives

*   Define "aggregation" and give examples of its use.
*   Write queries that compute aggregated values.
*   Trace the execution of a query that performs aggregation.
*   Explain how missing data is handled during aggregation.
</div>


We now want to calculate ranges and averages for our data.
We know how to select all of the dates from the `Visited` table:


<pre class="in"><code>%load_ext sqlitemagic</code></pre>


<pre class="in"><code>%%sqlite survey.db
select dated from Visited;</code></pre>

<div class="out"><table>
<tr><td>1927-02-08</td></tr>
<tr><td>1927-02-10</td></tr>
<tr><td>1939-01-07</td></tr>
<tr><td>1930-01-12</td></tr>
<tr><td>1930-02-26</td></tr>
<tr><td>None</td></tr>
<tr><td>1932-01-14</td></tr>
<tr><td>1932-03-22</td></tr>
</table></div>


but to combine them,
wee must use an [aggregation function](../../gloss.html#aggregation-function)
such as `min` or `max`.
Each of these functions takes a set of records as input,
and produces a single record as output:


<pre class="in"><code>%%sqlite survey.db
select min(dated) from Visited;</code></pre>

<div class="out"><table>
<tr><td>1927-02-08</td></tr>
</table></div>


<img src="img/sql-aggregation.svg" alt="SQL Aggregation" />


<pre class="in"><code>%%sqlite survey.db
select max(dated) from Visited;</code></pre>

<div class="out"><table>
<tr><td>1939-01-07</td></tr>
</table></div>


`min` and `max` are just two of
the aggregation functions built into SQL.
Three others are `avg`,
`count`,
and `sum`:


<pre class="in"><code>%%sqlite survey.db
select avg(reading) from Survey where quant=&#39;sal&#39;;</code></pre>

<div class="out"><table>
<tr><td>7.20333333333</td></tr>
</table></div>


<pre class="in"><code>%%sqlite survey.db
select count(reading) from Survey where quant=&#39;sal&#39;;</code></pre>

<div class="out"><table>
<tr><td>9</td></tr>
</table></div>


<pre class="in"><code>%%sqlite survey.db
select sum(reading) from Survey where quant=&#39;sal&#39;;</code></pre>

<div class="out"><table>
<tr><td>64.83</td></tr>
</table></div>


We used `count(reading)` here,
but we could just as easily have counted `quant`
or any other field in the table,
or even used `count(*)`,
since the function doesn't care about the values themselves,
just how many values there are.

SQL lets us do several aggregations at once.
We can,
for example,
find the range of sensible salinity measurements:


<pre class="in"><code>%%sqlite survey.db
select min(reading), max(reading) from Survey where quant=&#39;sal&#39; and reading&lt;=1.0;</code></pre>

<div class="out"><table>
<tr><td>0.05</td><td>0.21</td></tr>
</table></div>


We can also combine aggregated results with raw results,
although the output might surprise you:


<pre class="in"><code>%%sqlite survey.db
select person, count(*) from Survey where quant=&#39;sal&#39; and reading&lt;=1.0;</code></pre>

<div class="out"><table>
<tr><td>lake</td><td>7</td></tr>
</table></div>


Why does Lake's name appear rather than Roerich's or Dyer's?
The answer is that when it has to aggregate a field,
but isn't told how to,
the database manager chooses an actual value from the input set.
It might use the first one processed,
the last one,
or something else entirely.

Another important fact is that when there are no values to aggregate,
aggregation's result is "don't know"
rather than zero or some other arbitrary value:


<pre class="in"><code>%%sqlite survey.db
select person, max(reading), sum(reading) from Survey where quant=&#39;missing&#39;;</code></pre>

<div class="out"><table>
<tr><td>None</td><td>None</td><td>None</td></tr>
</table></div>


One final important feature of aggregation functions is that
they are inconsistent with the rest of SQL in a very useful way.
If we add two values,
and one of them is null,
the result is null.
By extension,
if we use `sum` to add all the values in a set,
and any of those values are null,
the result should also be null.
It's much more useful,
though,
for aggregation functions to ignore null values
and only combine those that are non-null.
This behavior lets us write our queries as:


<pre class="in"><code>%%sqlite survey.db
select min(dated) from Visited;</code></pre>

<div class="out"><table>
<tr><td>1927-02-08</td></tr>
</table></div>


instead of always having to filter explicitly:


<pre class="in"><code>%%sqlite survey.db
select min(dated) from Visited where dated is not null;</code></pre>

<div class="out"><table>
<tr><td>1927-02-08</td></tr>
</table></div>


Aggregating all records at once doesn't always make sense.
For example,
suppose Gina suspects that there is a systematic bias in her data,
and that some scientists' radiation readings are higher than others.
We know that this doesn't work:


<pre class="in"><code>%%sqlite survey.db
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant=&#39;rad&#39;;</code></pre>

<div class="out"><table>
<tr><td>roe</td><td>8</td><td>6.56</td></tr>
</table></div>


because the database manager selects a single arbitrary scientist's name
rather than aggregating separately for each scientist.
Since there are only five scientists,
she could write five queries of the form:


<pre class="in"><code>%%sqlite survey.db
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant=&#39;rad&#39;
and   person=&#39;dyer&#39;;</code></pre>

<div class="out"><table>
<tr><td>dyer</td><td>2</td><td>8.81</td></tr>
</table></div>


but this would be tedious,
and if she ever had a data set with fifty or five hundred scientists,
the chances of her getting all of those queries right is small.

What we need to do is
tell the database manager to aggregate the hours for each scientist separately
using a `group by` clause:


<pre class="in"><code>%%sqlite survey.db
select   person, count(reading), round(avg(reading), 2)
from     Survey
where    quant=&#39;rad&#39;
group by person;</code></pre>

<div class="out"><table>
<tr><td>dyer</td><td>2</td><td>8.81</td></tr>
<tr><td>lake</td><td>2</td><td>1.82</td></tr>
<tr><td>pb</td><td>3</td><td>6.66</td></tr>
<tr><td>roe</td><td>1</td><td>11.25</td></tr>
</table></div>


`group by` does exactly what its name implies:
groups all the records with the same value for the specified field together
so that aggregation can process each batch separately.
Since all the records in each batch have the same value for `person`,
it no longer matters that the database manager
is picking an arbitrary one to display
alongside the aggregated `reading` values.


Just as we can sort by multiple criteria at once,
we can also group by multiple criteria.
To get the average reading by scientist and quantity measured,
for example,
we just add another field to the `group by` clause:


<pre class="in"><code>%%sqlite survey.db
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
group by person, quant;</code></pre>

<div class="out"><table>
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
</table></div>


Note that we have added `person` to the list of fields displayed,
since the results wouldn't make much sense otherwise.

Let's go one step further and remove all the entries
where we don't know who took the measurement:


<pre class="in"><code>%%sqlite survey.db
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
where    person is not null
group by person, quant
order by person, quant;</code></pre>

<div class="out"><table>
<tr><td>dyer</td><td>rad</td><td>2</td><td>8.81</td></tr>
<tr><td>dyer</td><td>sal</td><td>2</td><td>0.11</td></tr>
<tr><td>lake</td><td>rad</td><td>2</td><td>1.82</td></tr>
<tr><td>lake</td><td>sal</td><td>4</td><td>0.11</td></tr>
<tr><td>lake</td><td>temp</td><td>1</td><td>-16.0</td></tr>
<tr><td>pb</td><td>rad</td><td>3</td><td>6.66</td></tr>
<tr><td>pb</td><td>temp</td><td>2</td><td>-20.0</td></tr>
<tr><td>roe</td><td>rad</td><td>1</td><td>11.25</td></tr>
<tr><td>roe</td><td>sal</td><td>2</td><td>32.05</td></tr>
</table></div>


Looking more closely,
this query:

1.  selected records from the `Survey` table
    where the `person` field was not null;

2.  grouped those records into subsets
    so that the `person` and `quant` values in each subset
    were the same;

3.  ordered those subsets first by `person`,
    and then within each sub-group by `quant`;
    and

4.  counted the number of records in each subset,
    calculated the average `reading` in each,
    and chose a `person` and `quant` value from each
    (it doesn't matter which ones,
    since they're all equal).


#### Challenges

1.  How many temperature readings did Frank Pabodie record,
    and what was their average value?

2.  The average of a set of values is the sum of the values
    divided by the number of values.
    Does this mean that the `avg` function returns 2.0 or 3.0
    when given the values 1.0, `null`, and 5.0?

3.  We want to calculate the difference between
    each individual radiation reading
    and the average of all the radiation readings.
    We write the query:

    ~~~
    select reading - avg(reading) from Survey where quant='rad';
    ~~~

    What does this actually produce, and why?

4.  The function `group_concat(field, separator)`
    concatenates all the values in a field
    using the specified separator character
    (or ',' if the separator isn't specified).
    Use this to produce a one-line list of scientists' names,
    such as:

    ~~~
    William Dyer, Frank Pabodie, Anderson Lake, Valentina Roerich, Frank Danforth
    ~~~

    Can you find a way to order the list by surname?


<div class="keypoints" markdown="1">
#### Key Points

*   An aggregation function combines many values to produce a single new value.
*   Aggregation functions ignore `null` values.
*   Aggregation happens after filtering.
</div>
