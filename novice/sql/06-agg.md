
## Aggregation

#### Objectives

*   Define "aggregation" and give examples of its use.
*   Write queries that compute aggregated values.
*   Trace the execution of a query that performs aggregation.
*   Explain how missing data is handled during aggregation.

We now want to calculate ranges and averages for our data.
We know how to select all of the dates from the `Visited` table:

In[1]:

```
%load_ext sqlitemagic
```

In[2]:

```
%%sqlite survey.db
select dated from Visited;
```




but to combine them,
wee must use an [aggregation function](../../gloss.html#aggregation-function)
such as `min` or `max`.
Each of these functions takes a set of records as input,
and produces a single record as output:

In[3]:

```
%%sqlite survey.db
select min(dated) from Visited;
```




<img src="files/img/sql-aggregation.svg" alt="SQL Aggregation" />

In[4]:

```
%%sqlite survey.db
select max(dated) from Visited;
```




`min` and `max` are just two of
the aggregation functions built into SQL.
Three others are `avg`,
`sum`,
and `count`:

In[5]:

```
%%sqlite survey.db
select avg(reading) from Survey where quant='sal';
```




In[6]:

```
%%sqlite survey.db
select sum(reading) from Survey where quant='sal';
```




In[7]:

```
%%sqlite survey.db
select count(reading) from Survey where quant='sal';
```




We used `count(reading)` here,
but since the function doesn't care about the values themselves,
just how many values there are, we could just as easily have
used count(`quant`). In practice, it's common
to use `count(*)` to count all rows.


SQL lets us do several aggregations at once.
We can,
for example,
find the range of sensible salinity measurements:

In[8]:

```
%%sqlite survey.db
select min(reading), max(reading) from Survey where quant='sal' and reading<=1.0;
```




Another important fact is that when there are no values to aggregate,
aggregation's result is "don't know"
rather than zero or some other arbitrary value. (Note that count returns a
reasonable value.)

In[9]:

```
%%sqlite survey.db
select count(reading), max(reading), sum(reading), avg(reading) from Survey where quant='missing';
```




One final important feature of aggregation functions is that
they ignore nulls. Note the full results below followed by the aggregated
counts.

In[10]:

```
%%sqlite survey.db
select ident, dated from Visited;
```




In[11]:

```
%%sqlite survey.db
select count(ident), count(dated) from Visited;
```




In constrast to non-aggregate functions/operations
(sometimes called scalar functions) where any null in the input results in a
null,
aggregate functions simply ignore nulls. While inconsistent with scalar
functions,
this behavior actually makes aggregate functions simpler and more useful.

The select list of our queries above contained only aggergate functions.
SQLite allows us to combine aggregated results with raw results,
but *be careful*: the results can be suprising.

In[12]:

```
%%sqlite survey.db
select person, count(*) from Survey where quant='sal' and reading<=1.0;
```




Is this what you expected?

The aggregate function `count(*)` in the select list cues the database manager
that you want to summarize all the results into a single row; but you also asked
for a `person`. Which person should it choose to represent all the original
rows?
In this example, why does Lake's name appear rather than Roerich's or Dyer's?
The answer is that different database managers handle this differently: SQLite
arbitrarily chooses a value from the original results; some database managers
present a null, and some database managers return an error. Fortunately, there
are several techniques to meaningfully combine aggregate results and column
values. Consider this example:

*Gina suspects that there is a systematic bias in her data, and that some
scientists' radiation readings are higher than others.
She would like to see the counts and averages for each scientist individually.*

From above, we know that this won't answer Gina's question:

In[13]:

```
%%sqlite survey.db
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant='rad';
```




because the database manager selects a single arbitrary scientist's name
rather than aggregating separately for each scientist.
Since there are only five scientists,
she could write five queries of the form:

In[14]:

```
%%sqlite survey.db
select person, count(reading), round(avg(reading), 2)
from  Survey
where quant='rad'
and   person='dyer';
```




but this would be tedious,
and if she ever had a data set with fifty or five hundred scientists,
the chances of her getting all of those queries right is small.

What we need to do is
tell the database manager to aggregate the hours for each scientist separately
using a `group by` clause:

In[15]:

```
%%sqlite survey.db
select   person, count(reading), round(avg(reading), 2)
from     Survey
where    quant='rad'
group by person;
```




`group by` does exactly what its name implies:
groups all the records with the same value for the specified field together
so that aggregation can process each batch separately.
Since all the records in each batch have the same value for `person`,
the database manager uses that value alongside the aggregated `reading` values.

Just as we can sort by multiple criteria at once,
we can also group by multiple criteria.
To get the average reading by scientist and quantity measured,
for example,
we just add another field to the `group by` clause:

In[16]:

```
%%sqlite survey.db
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
group by person, quant;
```




Note that we have added `person` to the list of fields displayed,
since the results wouldn't make much sense otherwise.

Let's go one step further and remove all the entries
where we don't know who took the measurement:

In[17]:

```
%%sqlite survey.db
select   person, quant, count(reading), round(avg(reading), 2)
from     Survey
where    person is not null
group by person, quant
order by person, quant;
```




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

2. The query:
       select count(1), count(name) from Survey
    returned the row:
       21 | 19
    How many rows are in the Survey table?
    Can you explain why the counts don't match?

3.  The average of a set of values is the sum of the values
    divided by the number of values.
    Does this mean that the `avg` function returns 2.0 or 3.0
    when given the values 1.0, `null`, and 5.0?

4.  We want to calculate the difference between
    each individual radiation reading
    and the average of all the radiation readings.
    We write the query:

    ~~~
    select reading - avg(reading) from Survey where quant='rad';
    ~~~

    What does this actually produce, and why?

5.  The function `group_concat(field, separator)`
    concatenates all the values in a field
    using the specified separator character
    (or ',' if the separator isn't specified).
    Use this to produce a one-line list of scientists' names,
    such as:

    ~~~
    William Dyer, Frank Pabodie, Anderson Lake, Valentina Roerich, Frank
Danforth
    ~~~

    Can you find a way to order the list by surname?

#### Key Points

*   An aggregation function combines many values to produce a single new value.
*   Aggregation functions ignore `null` values.
*   Aggregation happens after filtering.
