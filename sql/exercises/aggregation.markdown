## Getting Started

In these exercises we will be using a database of that contains information
about Nobel Prizes:  [nobel_prizes.sqlite](nobel_prizes.sqlite). 

## The exercises

**Question:** Fetch the earliest year recorded in the database.

**Answer:**
> ```sql
> SELECT MIN(Year) FROM Nobel_Prizes;
> ```

---

**Question:** Fetch the number of records in the database.

**Answer:**
> ```sql
> SELECT COUNT(*) FROM Nobel_Prizes;
> ```
	
---

**Question:** Fetch the area, year, and number of awards given out for each area
in each year.

**Answer:**
> ```sql
> SELECT Area, Year, COUNT(*) FROM Nobel_Prizes GROUP BY Area, Year;
> ```
	
---

**Question:** Fetch the area, and number of awards given out in each area,
sorted in ascending order by area.
   
**Answer:**
> ```sql
> SELECT Area, COUNT(*) 
> FROM Nobel_Prizes 
> GROUP BY Area 
> ORDER BY Area ASC;
> ```
	
---

**Question:** Fetch the area, and number of awards given out in each area,
sorted in descending order by number of awards.

**Answer:**
> 
> ```sql
> SELECT Area, COUNT(*) 
> FROM Nobel_Prizes 
> GROUP BY Area 
> ORDER BY Count(*) DESC;
> ```

---

**Question:** Fetch the area, and number of awards given out in each area,
except for in chemistry.

**Answer:**
> 
> ```sql
> SELECT Area, COUNT(*) 
> FROM Nobel_Prizes 
> WHERE Area != "chemistry"
> GROUP BY Area;
> ```
