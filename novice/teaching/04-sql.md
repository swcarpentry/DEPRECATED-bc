---
layout: lesson
root: ../..
title: Using Databases and SQL
level: novice
---
Relational databases are not as widely used in science as in business,
but they are still a common way to store large data sets with complex structure.
Even when the data itself isn't in a database,
the metadata could be:
for example,
meteorological data might be stored in files on disk,
but data about when and where observations were made,
data ranges,
and so on could be in a database
to make it easier for scientists to find what they want to.

#### Teaching Notes

*   The first few sections (up to "Missing Data") usually go very quickly.
    The pace usually slows down a bit when null values are discussed
    mostly because learners have a lot of details to keep straight by this point.
    Things *really* slow down during the discussion of joins,
    but this is the key idea in the whole lesson:
    important ideas like primary keys and referential integrity
    only make sense once learners have seen how they're used in joins.
    It's worth going over things a couple of times if necessary (with lots of examples).

*   The sections on creating and modifying data,
    and programming with databases,
    can be dropped if time is short.
    Of the two,
    people seem to care most about how to add data (which only takes a few minutes to demonstrate).


*   Overall,
    this material takes three hours to present assuming that a short exercise is done with each topic.

*   Simple calculations are actually easier to do in a spreadsheet;
    the advantages of using a database become clear as soon as filtering and joins are needed.
    Instructors may therefore want to show a spreadsheet with the information from the four database tables
    consolidated into a single sheet,
    and demonstrate what's needed in both systems to answer questions like,
    "What was the average radiation reading in 1931?"

*   Some learners may have heard that NoSQL databases
    (i.e., ones that don't use the relational model)
    are the next big thing,
    and ask why we're not teaching those.
    The answers are:
    1.  Relational databases are far more widely used than NoSQL databases.
    2.  We have far more experience with relational databases than with any other kind,
        so we have a better idea of what to teach and how to teach it.
    3.  NoSQL databases are as different from each other as they are from relational databases.
        Until a leader emerges, it isn't clear *which* NoSQL database we should teach.

*   Run `sqlite3 survey.db < gen-survey-database.sql`
    to re-create survey database before loading notebooks.
