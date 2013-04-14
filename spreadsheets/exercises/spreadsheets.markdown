## Exercise 1

Use the spreadsheet [spreadsheets_ex1.xls](spreadsheets_ex1.xls) for this exercise.

This spreadsheet contains a 3 x 3 matrix of numbers from one to nine.

1. In the orange cell, write a formula that adds up all the numbers above it. 
   **Answer:**
   `=SUM(B2:B4)`

2. Copy the orange cell and paste it into the green cell.  What value did it
produce?  Why?
   **Answer:**
   > The green cell should contain '15'.  When you copy a formula that uses
   > relative references, each reference is shifted over when you copy and
   > paste.  If you paste the formula one column to the left of the original,
   > all the references are translated one column to the left.

3. Copy the orange cell and paste it into the blue cell.  What value did it
produce?  Why?
   **Answer:**
   > The value should be '9'.  Since the blue cell is two cells to the left and
   > two cells below the orange, the references are translated the same
   > distance.  The formula ends up adding up cells D4 to D6, which is '9' and
   > two empty cells.


## Exercise 2

Use the spreadsheet [spreadsheets_ex2.xls](spreadsheets_ex2.xls) for this exercise.

In the orange cell, write a formula that applies the tax in the blue cell to the
cost in the pink cell.  Design the formula so that when you select the orange
and green cells and use the Fill -> Down command, the appropriate tax is
computed for each item.

   **Answer:**
   `=$E$2 * B2`


## Exercise 3

Use the spreadsheet [spreadsheets_ex3.xls](spreadsheets_ex3.xls) for this exercise.

Safety is a growing concern in your lab.  Employees have claimed that workplace
injuries have increased ever since the extra large quantum bogon collider was
turned on 4 months ago.  You would like to dispute these claims by showing that
on a monthly basis, injuries have been decreasing.  The spreadsheet below lists
injuries by department.  Create a pivot table which gives a breakdown of the
total number of injuries by month.

   **Hint:**
   > You can drag any of the fields into the "Drop Data Items Here"

   **Answer:**
   > After creating an empty pivot table, drag the "Month" field into either the
   > "Drop Column Fields Here" or "Drop Row Fields Here" areas.  Then drag any
   > column into "Drop Data Items Here".  If you pick any of the numeric fields,
   > you'll need change the field settings so that it summarizes by "Count"
   > instead of "Sum".


## Exercise 4

Create a column chart of the number of injuries by month from Exercise 3.
Months should go along the X-axis and injuries along Y-axis.

**Answer:**
> Select the cells containing the month names and the totals from the pivot
> table in Exercise 3.  Next, create a new clustered column chart from your
> selection.
