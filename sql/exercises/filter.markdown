## Getting Started

In these exercises we will be using a database of that contains information
about Nobel Prizes: [nobel_prizes.sqlite](nobel_prizes.sqlite). 

## The exercises

---

**Question:** Find all of the peace prizes awarded before 1910.

**Answer:**
> ```sql 
> SELECT * FROM Nobel_Prizes WHERE Year < 1910 AND Area = 'peace';
> ```
	
---

**Question:** Find all of the prizes awarded to Marie Curie:

**Answer:**
> ```sql
> SELECT * FROM Nobel_Prizes WHERE Name = 'Marie Curie';
> ```
> 
> Notice that we put quotation marks around Marie Curie to indicate that it is a
> string.
	
---

**Question:** Find all of the prizes that were awarded in areas of chemistry, peace, and
   physics only.  

**Answer:**
> Here are two ways to do this.  The first way uses the OR clause to match
> either 'chemistry', 'peace', or 'physics' in the Area field:
>
> ```sql
> SELECT * FROM Nobel_Prizes
> WHERE Area = 'chemistry' OR Area = 'peace' OR Area = 'physics';
> ```
>
> The other way to write this is by using an IN clause:
>
> ```sql
> SELECT * FROM Nobel_Prizes
> WHERE Area IN ('chemistry', 'peace', 'physics');
> ```
>
>  Both of these queries are equivalent and return the same results.
	
---

**Question:** Find all of the prizes that were awarded in areas _other than_ Chemistry,
   Peace, and Physics.  

**Answer:**
> Again, there are two ways to do this.  The first way uses the AND clause to
> say that neither 'chemistry', 'peace', nor 'physics' should occur the Area
> field:
> 
> ```sql
> SELECT * FROM Nobel_Prizes
> WHERE Area != 'chemistry' AND Area != 'peace' AND Area != 'physics';
> ```
>
> The other way to write this is by using NOT IN clause:
>
> ```sql
> SELECT * FROM Nobel_Prizes
> WHERE Area NOT IN ('chemistry', 'peace', 'physics');
> ```
>
> Again, both of these queries are equivalent but using the NOT IN clause makes
> the query much more succinct.
