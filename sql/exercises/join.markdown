## Getting Started

In these exercises we will be using a database of that contains information
about experiments: [experiments.sqlite](experiments.sqlite). 

## The exercises

In these exercises we use a different database than the previous exercises.  The
database contains three tables which track data about the projects and their
experiments being conducted in a laboratory, and the scientists involved.
	
- **`Person`** -- has an entry for each scientist in the lab.
- **`Experiment`** -- has an entry for each project and experiment id.
- **`ExperimentDetails`** -- has an entry for each project, and experiment and
  lists other details about the project itself.

---

**Question:** Display all of the information from both the `Experiment `and
`ExperimentDetail `table for each experiment.Â  That is, match rows where the
"Project" field is the same.

**Answer:**
> ```sql
> SELECT *
> FROM Experiment
> JOIN ExperimentDetail
> ON Experiment.Project = ExperimentDetail.Project;
> ```

---

**Question:** Again, join the `Experiment `table with the `ExperimentDetail
`table, returning only the rows where "Project" field is the same.  This time
though, use the alias "e" for the `Experiment `table, and "ed" for the
`ExperimentDetail` table.

**Answer:**
> ```sql
> SELECT *
> FROM Experiment e
> JOIN ExperimentDetail ed
> ON e.Project = ed.Project;
> ```

---

**Question:** Again, using aliases, join the `Experiment `table with the
`ExperimentDetail `table, returning only the rows where "Project" field is the
same.  This time, only fetch the "LoginID" and "Project" fields from the
`Experiment `table, and the "ExperimentName" field from `ExperimentDetail`.

**Answer:**
> ```sql
> SELECT e.LoginID, e.Project, ed.ExperimentName
> FROM Experiment e
> JOIN ExperimentDetail ed
> ON e.Project = ed.Project;
> ```

---

**Question:** Join the `Person`, `Experiment`, and `ExperimentDetail `tables.
Join the `Person `and `Experiment `tables on the "LoginID" field, and the
`Experiment `and `ExperimentDetail `tables on the "Project" AND "Experiment"
fields.  Fetch all the fields, and use aliases for the tables.

**Answer:**
> ```sql
> SELECT *
> FROM Person p
> JOIN Experiment e
> ON p.LoginID = e.loginID
> JOIN ExperimentDetail ed
> ON (e.Project = ed.Project AND e.Experiment = ed.Experiment);
> ```
