---
layout: lesson
root: ../..
---

## Missing Data


<div class="objectives" markdown="1">
#### Objectives

*   Explain how databases represent missing information.
*   Explain the three-valued logic databases use when manipulating missing information.
*   Write queries that handle missing information correctly.
</div>


Real-world data is never complete&mdash;there are always holes.
Databases represent these holes using special value called `null`.
`null` is not zero, `False`, or the empty string;
it is a one-of-a-kind value that means "nothing here".
Dealing with `null` requires a few special tricks
and some careful thinking.

To start,
let's have a look at the `Visited` table.
There are eight records,
but #752 doesn't have a date&mdash;or rather,
its date is null:


<pre class="in"><code>%load_ext sqlitemagic</code></pre>


<pre class="in"><code>%%sqlite survey.db
select * from Visited;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>752</td><td>DR-3</td><td>None</td></tr>
<tr><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></div>


Null doesn't behave like other values.
If we select the records that come before 1930:


<pre class="in"><code>%%sqlite survey.db
select * from Visited where dated&lt;&#39;1930-00-00&#39;;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
</table></div>


we get two results,
and if we select the ones that come during or after 1930:


<pre class="in"><code>%%sqlite survey.db
select * from Visited where dated&gt;=&#39;1930-00-00&#39;;</code></pre>

<div class="out"><table>
<tr><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></div>


we get five,
but record #752 isn't in either set of results.
The reason is that
`null<'1930-00-00'`
is neither true nor false:
null means, "We don't know,"
and if we don't know the value on the left side of a comparison,
we don't know whether the comparison is true or false.
Since databases represent "don't know" as null,
the value of `null<'1930-00-00'`
is actually `null`.
`null>='1930-00-00'` is also null
because we can't answer to that question either.
And since the only records kept by a `where`
are those for which the test is true,
record #752 isn't included in either set of results.

Comparisons aren't the only operations that behave this way with nulls.
`1+null` is `null`,
`5*null` is `null`,
`log(null)` is `null`,
and so on.
In particular,
comparing things to null with = and != produces null:


<pre class="in"><code>%%sqlite survey.db
select * from Visited where dated=NULL;</code></pre>

<div class="out"><table>

</table></div>


<pre class="in"><code>%%sqlite survey.db
select * from Visited where dated!=NULL;</code></pre>

<div class="out"><table>

</table></div>


To check whether a value is `null` or not,
we must use a special test `is null`:


<pre class="in"><code>%%sqlite survey.db
select * from Visited where dated is NULL;</code></pre>

<div class="out"><table>
<tr><td>752</td><td>DR-3</td><td>None</td></tr>
</table></div>


or its inverse `is not null`:


<pre class="in"><code>%%sqlite survey.db
select * from Visited where dated is not NULL;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>DR-1</td><td>1927-02-08</td></tr>
<tr><td>622</td><td>DR-1</td><td>1927-02-10</td></tr>
<tr><td>734</td><td>DR-3</td><td>1939-01-07</td></tr>
<tr><td>735</td><td>DR-3</td><td>1930-01-12</td></tr>
<tr><td>751</td><td>DR-3</td><td>1930-02-26</td></tr>
<tr><td>837</td><td>MSK-4</td><td>1932-01-14</td></tr>
<tr><td>844</td><td>DR-1</td><td>1932-03-22</td></tr>
</table></div>


Null values cause headaches wherever they appear.
For example,
suppose we want to find all the salinity measurements
that weren't taken by Dyer.
It's natural to write the query like this:


<pre class="in"><code>%%sqlite survey.db
select * from Survey where quant=&#39;sal&#39; and person!=&#39;lake&#39;;</code></pre>

<div class="out"><table>
<tr><td>619</td><td>dyer</td><td>sal</td><td>0.13</td></tr>
<tr><td>622</td><td>dyer</td><td>sal</td><td>0.09</td></tr>
<tr><td>752</td><td>roe</td><td>sal</td><td>41.6</td></tr>
<tr><td>837</td><td>roe</td><td>sal</td><td>22.5</td></tr>
</table></div>


but this query filters omits the records
where we don't know who took the measurement.
Once again,
the reason is that when `person` is `null`,
the `!=` comparison produces `null`,
so the record isn't kept in our results.
If we want to keep these records
we need to add an explicit check:


<pre class="in"><code>%%sqlite survey.db
select * from Survey where quant=&#39;sal&#39; and (person!=&#39;lake&#39; or person is null);</code></pre>

<div class="out"><table>
<tr><td>619</td><td>dyer</td><td>sal</td><td>0.13</td></tr>
<tr><td>622</td><td>dyer</td><td>sal</td><td>0.09</td></tr>
<tr><td>735</td><td>None</td><td>sal</td><td>0.06</td></tr>
<tr><td>752</td><td>roe</td><td>sal</td><td>41.6</td></tr>
<tr><td>837</td><td>roe</td><td>sal</td><td>22.5</td></tr>
</table></div>


We still have to decide whether this is the right thing to do or not.
If we want to be absolutely sure that
we aren't including any measurements by Lake in our results,
we need to exclude all the records for which we don't know who did the work.


#### Challenges

1.  Write a query that sorts the records in `Visited` by date,
    omitting entries for which the date is not known
    (i.e., is null).

1.  What do you expect the query:

    ~~~
    select * from Visited where dated in ('1927-02-08', null);
    ~~~

    to produce?
    What does it actually produce?

1.  Some database designers prefer to use
    a [sentinel value](../../gloss.html#sentinel-value)
    to mark missing data rather than `null`.
    For example,
    they will use the date "0000-00-00" to mark a missing date,
    or -1.0 to mark a missing salinity or radiation reading
    (since actual readings cannot be negative).
    What does this simplify?
    What burdens or risks does it introduce?


<div class="keypoints" markdown="1">
#### Key Points

*   Databases use `null` to represent missing information.
*   Any arithmetic or Boolean operation involving `null` produces `null` as a result.
*   The only operators that can safely be used with `null` are `is null` and `is not null`.
</div>
