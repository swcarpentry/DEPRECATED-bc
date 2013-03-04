## Overall Exercises


In the following exercises, suppose we have a database with the following two
tables:

- **Table Name:** Participants
- **Fields:** ParticipantID, Age, Gender, ConditionID, Score

- **Table Name:** Conditions
- **Fields:** ConditionID, Description

---

**Question:** Create the tables.  Run the following SQL statements (we haven't
covered SQL that modifies the structure of a database, but it should be
relatively clear how it works).
  
**Answer:**
> ```sql 
> CREATE TABLE Participants (
>    ParticipantID integer,
>    Age integer, 
>    Gender text, 
>    ConditionID integer, 
>    Score integer
> );
> 
> CREATE TABLE Condition (
>    ConditionID integer, 
>    Description text
> );
> ```

---

**Question:** Add data.  Run the following SQL statements.  You can also make up others.

**Answer:**
> ```sql
> INSERT INTO Participants VALUES (101, 22, 'm', 1, 56);
> INSERT INTO Participants VALUES (102, 20, 'f', 2, 113);
> INSERT INTO Participants VALUES (102, 20, 'f', 3, 95);
> INSERT INTO Participants VALUES (103, 23, 'm', 1, 87);
> INSERT INTO Participants VALUES (103, 23, 'm', 5, 90);
> INSERT INTO Participants VALUES (104, 22, 'm', 4, 92);
> INSERT INTO Participants VALUES (105, 28, 'f', 3, 3);
> INSERT INTO Participants VALUES (106, NULL, 'm', 4, 70);
> INSERT INTO Participants VALUES (106, NULL, 'm', 5, 100);
> 
> INSERT INTO Condition VALUES (1, 'Red lights only, silence.');
> INSERT INTO Condition VALUES (2, 'Red lights only, noise.');
> INSERT INTO Condition VALUES (3, 'Blue lights only, silence.');
> INSERT INTO Condition VALUES (4, 'Blue lights only, noise.');
> INSERT INTO Condition VALUES (5, 'Blue and red lights, silence.');
> INSERT INTO Condition VALUES (6, 'Blue and red lights, noise.');
> ```

---

**Question:** Write a query to fetch the "ParticipantID" and "Age" fields from
the Participants table. 

**Answer:**
> ```sql
> SELECT ParticipantID, Age FROM Participants;
> ```

---

**Question:** Write a query to fetch the distinct ages in the Participants table. 

**Answer:**
> ```sql
> SELECT DISTINCT Age FROM Participants;
> ```

---

**Question:** Suppose that "Score" field in the Participants table is the
participant's raw score out of 120.  Write a query to fetch the ParticipantIDs
and scores, but with the scores expressed out of 100 instead of 120, and rounded
to one decimal place. 

**Answer:**
> ```sql
> SELECT ParticipantID, ROUND(Score / 120.0 * 100, 1) FROM Participants;
> ```

---

**Question:** Write a query to fetch the ParticipantID and Age of all female
participants who scored higher than 50.

**Answer:**
> ```sql
> SELECT ParticipantID, Age
> FROM Participants
> WHERE Gender = "f" AND Score > 50;
> ```

---

**Question:** Write a query to fetch the ParticipantID of all participants who
were not in condition 1, 3, or 5. 

> ```sql
> SELECT ParticipantID
> FROM Participants
> WHERE ConditionID NOT IN (1, 3, 5);
> ```

---

**Question:** Write a query to fetch all the records in the Participants table,
sorted in ascending order by age.

**Answer:**
> ```sql
> SELECT * FROM Participants ORDER BY Age ASC;
> ```

---

**Question:** Write a query to fetch all the records in the Participants table,
sorted randomly. 

**Answer:**
> ```sql
> SELECT * FROM Participants ORDER BY RANDOM();
> ```

---

**Question:** Write a query to get the average score for each gender. 

**Answer:**
> ```sql
> SELECT Gender, AVG(Score) FROM Participants GROUP BY Gender;
> ```

---

**Question:** Write a query to join the Participants and Conditions tables on
the field "ConditionID".  Fetch all the fields, and use aliases for your tables. 

**Answer:**
> ```sql
> SELECT *
> FROM Participants p
> JOIN Condition c
> ON p.ConditionID = c.ConditionID;
> ```

---

**Question:** Suppose we want to fetch all the records in the Participants table
for which the "Gender" field is missing data.  Why will the following query not
work, and what do we have to change to make it work?

**Answer:**
> ```sql
> SELECT * FROM Participants WHERE Gender = NULL; 
> ```
>
> NULL is a special value that cannot be compared using the usual operators.
> We need to use "IS" instead of "=" to make the query work.

--- 

**Question:** Write a query to count the number of students for which the gender
*IS* recorded. 

**Answer:**
> ```sql
> SELECT COUNT(*) FROM Participants WHERE Gender IS NOT NULL;
> ```

---

**Question:** Suppose that the experiment has used a within-subjects design,
such that participants may have been in more than one condition.  Write a query
to fetch the IDs of the participants who did not take part in condition 3. 

**Answer:**
> ```sql
> SELECT DISTINCT ParticipantID FROM Participants
> WHERE ParticipantID NOT IN
>   (SELECT DISTINCT ParticipantID FROM Participants
>    WHERE ConditionID = 3);
> ```
