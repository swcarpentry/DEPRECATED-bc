## Getting Started


In these exercises we will be using a database of that contains information
about Nobel Prizes: [nobel_prizes.sqlite](nobel_prizes.sqlite). 

## The exercises
	
Question: Fetch all of the prizes, and sort them by Year.   

Answer:
> ```sql
> SELECT * FROM Nobel_Prizes ORDER BY Year;
> ```
> 
> You may have chosen to specify ASC or DESC after the Year column as well, like so:
> 
> ```sql
> SELECT * FROM Nobel_Prizes ORDER BY Year DESC;
> ```

---

Question: Fetch all of the prizes, and sort them by Area in descending order,
and within each area sort by Year in ascending order.  

Answer: 
> ```sql
> SELECT * FROM Nobel_Prizes   
> ORDER BY Area DESC, Year ASC;
> ```

---

Question: Fetch all of the Nobel Prizes awarded in or before 1950 and sort the
results by the Area.     

Answer: 
> ```sql
> SELECT * FROM Nobel_Prizes
> WHERE Year <= 1950
> ORDER BY Area;
> ```
> 
> Notice that we put the ORDER BY clause _after_ the WHERE clause. In SQL, the
> different clauses have to appear in a certain order.  We'll discuss this again
> in a future lecture.
