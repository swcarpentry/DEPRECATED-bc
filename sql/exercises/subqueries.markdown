## Getting Started

In these exercises we will be using a database of that contains information
about experiments:  [experiments_subqueries.sqlite](experiments_subqueries.sqlite). 


## The exercises

---

**Question:** From the `Experiment` table, fetch the LoginIDs of scientists who
have **not** worked on the Teleportation project.

One approach is to see the solution as the being to gather all login IDs, and
then filter out those that represent scientists who have worked on the
Teleportation project.

Why doesn't this work? 

```sql
SELECT DISTINCT LoginID
FROM Experiment
WHERE Project != 'Teleportation';
```

**Hint:** 
> Consider following these steps in order to write the query:
> 
> 1. Write a query that fetches all of the LoginIDs (without duplicates)
> 2. Write a query that fetches all of the LoginIDs representing scientists who
>    have worked on the Teleportation project.
> 3. Combine the two queries so that you filter out those IDs returned from the
>    second query from the first.


**Answer:** 
> ```sql
> SELECT DISTINCT LoginID FROM Experiment
> WHERE LoginID NOT IN
> (SELECT DISTINCT LoginID FROM Experiment
> WHERE Project = 'Teleportation');
> ```sql

---

**Question:** From the`Experiment` table, fetch the LoginIDs of scientists who
have **not** worked on any experiments that took more than six hours.


**Hint:** 
> You don't need a hint!  This is very similar to the previous exercise.


**Answer:** 
> ```sql
> SELECT DISTINCT LoginID FROM Experiment 
> WHERE LoginID NOT IN 
>   (SELECT DISTINCT LoginID FROM Experiment 
>    WHERE Hours > 6); 
> ```sql

---

**Question:** From the `Experiment` table, fetch the LoginIDs of the scientists
who have done experiments on at least two different dates.


**Hint:** 
> One approach to this problem is to think of the solution like this. For a
> record to be included in the results, there must be another record in the
> table with the same `LoginID`, but with a different `ExperimentDate`.
> 
> You can perform this query by writing an an outer query that fetches all of
> the LoginIDs, and a nested query that checks that a record exists with that
> same LoginID but with a different experiment date. 
 
**Answer:**
> ```sql
> SELECT DISTINCT LoginID FROM Experiment e1
> WHERE LoginID IN
>   (SELECT DISTINCT LoginID FROM Experiment e2
>    WHERE e1.LoginID = e2.LoginID AND e1.ExperimentDate != e2.ExperimentDate);
> ```sql

---

**Question:** From the `Experiment` table, fetch the LoginIDs of all the
scientists, and for each include the number of distinct dates on which they have
conducted experiments.

**Hint:** 
> One approach is to see the solution as a two stage process.  First, you get a
> list of all of the unique pairs of login IDs and experiment dates.  Then,
> given that list, you count up the experiment dates for each login ID. 
> 
> Start by just writing a query that fetches unique pairs of login ID and
> experiment dates.
> 
> ```sql
> SELECT DISTINCT LoginID, ExperimentDate FROM Experiment;
> ```

**Hint:** 
> You'll want to use a GROUP BY statement to aggregate the count over each login
> ID.

**Answer:** 
> ```sql
> SELECT LoginID, COUNT(*)
> FROM (SELECT DISTINCT LoginID, ExperimentDate FROM Experiment)
> GROUP BY LoginID;
> ```
