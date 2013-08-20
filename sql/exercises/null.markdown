## Getting Started


In these exercises we will be using a database of that contains information
about experiments: [experiments_null.sqlite](experiments_null.sqlite).

## The exercises

The Experiments database is described in the [JOIN
exercises](join.markdown). In the version of the database given above,
some of the data in the tables is missing, which you'll have to take into
account.

---

**Question:** Fetch the rows from the `Experiment` table that are missing
`ExperimentDate` data.

**Answer:**
> ```sql
> SELECT *
> FROM Experiment
> WHERE ExperimentDate IS NULL;
> ```
>
> You might have tried to write your WHERE-clause as "ExperimentDate = NULL".  Why
> doesn't that work?
>
> `NULL` is a special value.  Comparing it using the equality operator _always_
> returns false.
	
---

**Question:** Fetch the rows from the `Experiment` table that are **not**
missing `ExperimentDate` data.

**Answer:**
> ```sql
> SELECT *
> FROM Experiment
> WHERE ExperimentDate IS NOT NULL;
> ```

---

**Question:** Fetch the rows from the `Experiment` table that do not have the
value "1" in the `Experiment` field.

The rows you fetch must include rows where the `Experiment` field is some other
number other than 1.  But, as we've seen above, NULL is special so you'll
probably have to handle it explicitly.

**Answer:**
>
> ```sql
> SELECT * 
> FROM Experiment 
> WHERE Experiment != 1 OR Experiment IS NULL;
> ```
