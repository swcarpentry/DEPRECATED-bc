---
layout: lesson
root: ../..
title: Data mining text files on the unix command line
---

### Motivation
In many fields, you will encounter data stored in text files.  For complex text parsing and re-formatting tasks or well-known formats, it might be a good idea to use or write a shell script, python script, etc. that simplifies the job.  However, if you want to do a quick, ad-hoc data extraction, it can sometimes be easier to use some simple, command-line  data mining techniques to get what you need in less time than it would take to write a script.


### Objective

Learn to use a simple set of unix command line tools for "quick and dirty" data extraction from two commonly encountered text file types.

### Level

This is a novice-to-intermediate level lesson.  It is assumed that you know the unix command line basics and have seen a few of the commands before.

### Your toolkit

This is a list of a few commands that we will use in the examples below. 

<table>
  <tr><th>|</th><td><i>strings together the inputs/outputs of a series of commands</i></td></tr>
  <tr><th>grep</th><td><i>searches for patterns in text</i></td></tr>
  <tr><th>sort</th><td><i>orders lines in text</i></td></tr>
  <tr><th>head</th><td><i>prints the top N lines of text</i></td></tr>
  <tr><th>tail</th><td><i>prints the bottom N lines of text</i></td></tr>
  <tr><th>uniq</th><td><i>reports or filters repeated lines of text</i></td></tr>
  <tr><th>wc</th><td><i>counts words, characters or lines</i></td></tr>
  <tr><th>bc</th><td><i>does simple math on the command line</i></td></tr>
</table>


### Example One: gene expression data

You have a tab-delimited text file, gene_exp.txt, that contains data from a differential gene expression analysis.  Each line describes a comparison of numerical expression levels for one gene in two samples.

#### What does the file look like 

What is the file structure? Without options, ***head*** will print the top 10 lines of the file

<pre>
$ head  gene_exp.txt
gene	sample_1	sample_2	status	value_1	value_2	significant
AT1G01010		WT		hy5	NOTEST	0	1.49367	no
AT1G01020		WT		hy5	NOTEST	7.27837	10.7195	no
AT1G01030		WT		hy5	NOTEST	1.18638	1.10483	no
AT1G01040		WT		hy5	NOTEST	0.239843	2.24208	no
AT1G01046		WT		hy5	NOTEST	0		0	no
AT1G01050		WT		hy5	OK	9.06975		23.5089	yes
AT1G01060		WT		hy5	NOTEST	4.04534		6.46964	no
AT1G01070		WT		hy5	NOTEST	1.24918		2.41377	no
AT1G01073		WT		hy5	NOTEST	0		0	no
</pre>

This is a well-formatted, tab-delimited text file.  The header line describes the columns for us.  We can use this to help answer some questions.  Note that ***sort*** and ***cut*** assume that the columns in each row are tab-delimited.

#### How many records are there in the file?

We can use ***wc -l*** to count the lines

<pre>
$ wc -l gene_exp.txt 
   33567 gene_exp.txt
</pre>


#### How many genes have enough data to perform the comparison (have 'OK' status)? How many had significantly different expression levels between samples?

We can search for *OK* and *yes* in the file, then count how many lines are returned by ***grep***.  Note the use of ***|***.  ***wc -l*** is acting on the text printed by ***grep***, not the input file.

<pre>
$ grep OK gene_exp.txt | wc -l
    4112
$ grep yes gene_exp.txt | wc -l
    1403
</pre>

For 33,567 genes, 4,112 had enough data to do a comparison and 1,403 had significantly different expression.


### Question:  What if the strings *yes* or *OK* appear in other columns of the file?

There is not a lot of room for free text in this example but it can't hurt to check.  This sort of thing happens all the time in real life!  We can use ***cut*** to remove the normal column (for example, column 7 for 'yes').  Remove column 7, then search for *yes*.

<pre>
$ cut -f1-6 gene_exp.txt | grep yes
</pre>

No results.  That is good.  For *OK* we need to search all columns except column 4:

<pre>
cut -f1-3,5-7 gene_exp.txt | grep OK
</pre>

No results.  The file is good.  Note that the ***-k*** argument could also have been expressed as ***-k1,2,3,5,6,7***


#### What are the 20 genes with the highest expression levels in sample 1 and differ significantly between samples?

We can use ***grep*** to get the 'yes' lines, then use ***sort*** to order the lines base on the numeric values in column 5 (value_1).  We pipe the output to ***head*** so we just look at the top 10 lines for now.  The ***k5*** argument means sort on column (key) 5; ***-n*** means sort numerically.  

<pre>
$ grep 'yes' gene_exp.txt | sort -k5 -n | head
AT1G40125    WT		  hy5	 OK  0	15.3962	yes
AT1G42040    WT		  hy5	 OK  0	23.5267	yes
AT1G42050    WT		  hy5	 OK  0	31.0539	yes
AT2G05915    WT		  hy5	 OK  0	61.649	yes
AT2G40802    WT		  hy5	 OK  0	551.414	yes
AT3G01345    WT		  hy5	 OK  0	29.1111	yes
AT3G22235    WT		  hy5	 OK  0	14.7018	yes
AT3G33073    WT		  hy5	 OK  0	18.2451	yes
AT3G42720    WT		  hy5	 OK  0	19.4265	yes
AT4G06530    WT		  hy5	 OK  0	20.3433	yes
</pre>

You may have noticed that cut and sort use different argumants for the same thing (column number).  The collection of tools in unix-like operating systems evolved over time from a variety of sources and authors, so their command line argumants are not always consistent.  If in doubt:

<pre>
man cut
</pre>

Note the the values in column 5 above are all zeros.  We are not quite there yet.  We can use the ***r*** flag to sort in descending order.

<pre>
$ grep 'yes' gene_exp.txt | sort -k5 -n -r | head
AT2G01021    WT		  hy5	 OK  282360  3.44931e+06	yes
AT1G08115    WT		  hy5	 OK  69434.3 35118.3		yes
ATCG00010    WT		  hy5	 OK  51851.7 23458.4		yes
ATCG00400    WT		  hy5	 OK  27712.3 2078		yes
AT5G41471    WT		  hy5	 OK  27289.6 3739.15		yes
XLOC_013786  WT		  hy5	 OK  24253.8 2509.04		yes
ATCG00630    WT		  hy5	 OK  22883.4 5164.64		yes
AT3G24615    WT		  hy5	 OK  13744.4 4184.12		yes
AT4G39363    WT		  hy5	 OK  11180.5 2136.42		yes
AT3G06895    WT		  hy5	 OK  9823.28 2143.76		yes
</pre>

OK, now we have the 10 most abundant genes in sample 1.  However, the question was "What are the 20 genes with the highest expression levels in sample 1...".  We can get the top 20 by adding the ***-20*** argument to ***head***.  Note also that we were asked for the gene, not the gene plus data.  We can use ***cut*** to extract what we want.  The ***-f1*** arument tells it to grab the first column.

<pre>
$ grep 'yes' gene_exp.txt | sort -k5 -n -r | head -20 | cut -f1
AT2G01021
AT1G08115
ATCG00010
ATCG00400
AT5G41471
XLOC_013786
ATCG00630
AT3G24615
AT4G39363
AT3G06895
XLOC_008330
ATCG00700
ATCG00390
AT3G41768
AT1G29930
AT1G79040
XLOC_001625
AT1G67090
AT3G56020
XLOC_032942
</pre>

And we have our answer!

### Example Two: quick and dirty web log analysis.

Suppose you are an admin on a high traffic web server and you are getting user complaints about the web site being very slow.  You check the database logs and find that there has been a spike in web traffic since about 11AM.  How can you check the apache access log to see what is going on?

Consider the log file access_log.  The exact format of the file can be customized in the web server configuration. The file contents are never pretty.

Your task is to find the top 10 users of the web site over the past hour.

#### How many lines in the file (ie, how many times was the web site accesses in the time period covered by this file)?

<pre>
$ wc -l access_log 
    37554 access_log
</pre>

#### What does the log file look like?

<pre>
$ head -5 access_log 
221.0.112.219 - - [07/Apr/2014:00:00:00 -0400] "GET /download/current/sql.gz HTTP/1.1" 200 63833565
58.95.175.121 - - [07/Apr/2014:00:00:11 -0400] "GET / HTTP/1.1" 200 22495
112.249.80.13 - - [07/Apr/2014:00:00:39 -0400] "HEAD /download/current/sql.gz HTTP/1.1" 200 -
165.246.204.254 - - [07/Apr/2014:00:01:23 -0400] "GET /ReactomeRESTfulAPI/RESTfulWS/queryById/Pathway/5388356 HTTP/1.1" 200 489
165.246.204.254 - - [07/Apr/2014:00:01:28 -0400] "POST /ReactomeRESTfulAPI/RESTfulWS/queryByIds/DatabaseObject HTTP/1.1" 200 490
</pre>

Not much to look at but we can see that the IP address of the browser is the first part of the record and that consistent time stamps are used.

#### It is 12 noon, April 7, 2014.  Isolate the records for the past hour.  

What time does the log file end?

<pre>
$ tail -1 access_log
adsl-4.46.190.51.tellas.gr - - [07/Apr/2014:11:59:57 -0400] "GET /cgi-bin/images/search.gif HTTP/1.1" 404 301
</pre>

It ends at noon.  We just need to know what line number should we start at.  Use ***grep*** to find the first line that start as 11AM.

<pre>
$ grep -n '07/Apr/2014:11:00' access_log | head -1
33580:66.249.74.215 - - [07/Apr/2014:11:00:02 -0400] "GET /img-tmp/650.4537338364293990.png HTTP/1.1" 404 308
</pre>

OK, the first line for the past hour is 33580.  How can we isolate only records from line 33580 and below?  We know the file has  37554 lines.  We can use ***bc*** and ***tail*** to get the lines we need.  For the first pass we will use ***head -1*** to make the the first line is the one we want.

<pre>
$ tail -3974 access_log | head -1 
146.185.30.53 - - [07/Apr/2014:11:00:03 -0400] "GET /cgi-bin/instancebrowser?DB=gk_current&ID=109581 HTTP/1.1" 200 68805
</pre>

Now we have the text isolated.  

#### Who is accessing the site and how often?  

Remember that the first section of each line is the the address of the browser.  The file is not tab-delimited, but we can use the ***-d' '*** argument for ***cut*** to tell it to use space as the delimiter.

Let's check first by piping the output to ***head***:

<pre>
$ tail -3974 access_log | cut -d' ' -f1 | head
146.185.30.53
5.10.83.18
90.83.115.106
90.83.115.106
90.83.115.106
90.83.115.106
90.83.115.106
90.83.115.106
90.83.115.106
90.83.115.106
</pre>

Looks good.  Now, how to we count the number of times each address occurs in the log?  ***uniq -c *** can report repeated lines.

<pre>
$ tail -3974 access_log | cut -d' ' -f1 |uniq -c |head 
   1 146.185.30.53
   1 5.10.83.18
   9 90.83.115.106
   1 mail.decharenton.fr
   4 90.83.115.106
   1 5.10.83.36
   1 5.10.83.40
  29 88.159.161.40
  16 14.98.38.62
   1 se-bravo.psycho.unibas.ch
</pre>

...but some of the addresses are listed more than once.  Let's try sorting before ***uniq***.

<pre>
$ tail -3974 access_log | cut -d' ' -f1 | sort | uniq -c | head 
   2 223.87.29.50
   1 217.69.133.232
  46 216.99.65.83
  24 213.1.212.91
   2 212.189.216.198
   1 211.136.10.53
   4 209.37.248.167
   3 208.115.113.86
  10 207.162.51.5
   1 206.83.48.110
</pre>


# Who are the the ten users between 11-12?
We can find out with another round of sorting at the end. Note that the last sort is numeric, in reverse (descending) order.

<pre>
$ tail -3974 access_log | cut -d' ' -f1 | sort | uniq -c | sort -n -r | head
 349 mrcbloxx.le.ac.uk
 270 se-bravo.psycho.unibas.ch
 256 hx-dnat-245.ebi.ac.uk
 190 193.63.220.131
 141 143.169.236.220
 129 77.125.74.117
 116 spelman-fw.spelman.edu
 105 91.103.43.50
 101 155.250.255.141
  99 hx-dnat-249.ebi.ac.uk
</pre>

Now that we know who is using the site most heavily, we can do some more directed searching, based on the addresses, to see what part of the web site they are hitting.
