---
layout: lesson
root: ../..
title: Summing Up
---

Every regular expressions is actually a little program: a piece of data
that tells a computer how to behave. The steps the computer goes through
when executing a Python or Fortran program are more complicated than the
ones it uses to compile and apply a regular expression, but the
principles are exactly the same.

One last point to take away from this chapter is that if we know we are
going to use regular expressions to read in data, we should choose a
format for that data that's easy for regular expressions to match.
Optional commas, tabs that might be repeated, and other things that make
data easy for people to type in actually make it harder for programs to
read that data reliably. This tension between what's easy for the
machine and what's easy for the user never goes away, but if we're
conscious of it, we can find a happy medium.
