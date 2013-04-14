## Getting Started

In these exercises we will be using a database of that contains information
about Nobel Prizes [nobel_prizes.sqlite](nobel_prizes.sqlite). 

This database only contains one table: `Nobel_Prizes`. 

## The exercises
	
---

**Question:** Fetch all of the records from the `Nobel_Prizes` table. 

> ```sql
> SELECT * FROM Nobel_Prizes;
> ```
> Notice how the `\*` is used to select all of the columns
> 
---

**Question:** Fetch all of the records from the `Nobel\_Prizes` table, but only
show return the award winner's name and year the prize was awarded. 

**Answer:**
> ```sql
> SELECT Year, Name FROM Nobel_Prizes;
> ```

---

**Question:** Fetch a list of all of the different areas the Nobel Prizes have
been awarded in. Don't return any duplicates! 

**Answer:**
> ```sql
> SELECT DISTINCT Area FROM Nobel_Prizes;
> ```

---

**Question:** The `UPPER(str)` function takes a string `str` and returns a new,
upper-case version of `str`. For instance, `SELECT UPPER("Hello");` returns the
string "HELLO".

Using `UPPER()`, write a query that fetches all of the Nobel Prize winner's
names, both as it is in the database, and in upper-case.

**Answer:**
> ```sql
> SELECT Name, UPPER(Name) FROM Nobel_Prizes;
> ```
