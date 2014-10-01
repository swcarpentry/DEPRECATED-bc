---
layout: lesson
root: ../..
---

# SQLite Setup 


This lesson will demonstrate how to install the example database to be used in the next chapters. To be able to follow the instructions, you need to know how to move around in your directory using the command line and how to execute commands from the command line. If you are not familiar with these topics, please read the [Unix Shell guide](http://software-carpentry.org/v5/novice/shell/index.html). In a later chapter, you will learn how to create and fill a database, but first we would like to demonstrate how an SQLite database works and therefore we provide the database up-front.


#### Objectives

*   Establish the example database to be used in the next chapters.
*   Check that the database is available and which tables are to be found.


## Installation


In order to execute the following lessons interactively, please install SQLite as mentioned in the [setup instructions](http://software-carpentry.org/v5/setup.html).

Then, generate a directory "software_carpentry_sql" at your chosen location, e.g.:


1) Open a command line terminal window.  
2) Type 


<pre class="in"><code>mkdir ~/swc/sql</code></pre>


3) Change into that directory 


<pre class="in"><code>cd ~/swc/sql</code></pre>


### How to download the "gen-survey-database.sql" file from github


After changing into the "~/swc/sql" directory, from within the directory you will now download the SQL file "gen-survey-database.sql" located on github: https://github.com/swcarpentry/bc/blob/master/novice/sql/gen-survey-database.sql  

Since this file is within a git repository, you will pull the single file locally without cloning the entire git repo. For this purpose, you can use either [GNU Wget](http://en.wikipedia.org/wiki/Wget), a command-line web-crawler software that supports downloading via HTTP, HTTPS, and FTP protocols; or use [cURL](http://en.wikipedia.org/wiki/CURL), a library and command-line tool for transferring data using various protocols. Both these tools are cross-platform and are supported for various operating systems. 

After installing Wget (or cURL) locally, run the following command on your terminal:

<i>[Tip: If you prefer to use cURL, replace "wget" with "curl -O" in the following command.]</i>


<pre class="in"><code>mom@durga:~/swc/sql$ wget https://raw.githubusercontent.com/swcarpentry/bc/master/novice/sql/gen-survey-database.sql</code></pre>


With the above command, Wget generates an HTTP request to pull the single raw file "gen-survey-database.sql" from github into your current directory and upon successful completion, your terminal will display the following output:


<pre class="in"><code>--2014-09-02 18:31:43--  https://raw.githubusercontent.com/swcarpentry/bc/master/novice/sql/gen-survey-database.sql
Resolving raw.githubusercontent.com (raw.githubusercontent.com)... 103.245.222.133
Connecting to raw.githubusercontent.com (raw.githubusercontent.com)|103.245.222.133|:443... connected.
HTTP request sent, awaiting response... 200 OK
Length: 3297 (3.2K) [text/plain]
Saving to: ‘gen-survey-database.sql’

100%[=========================================================================================================================&gt;] 3,297       --.-K/s   in 0.01s   

2014-09-02 18:31:45 (264 KB/s) - ‘gen-survey-database.sql’ saved [3297/3297]</code></pre>


Now that we have successfully pulled the single SQL file, lets generate the database "survey.db" and fill it according to the instructions with the data in "gen-survey-database.sql". To call the SQLite3 program, from the command line terminal, execute the following command:


<pre class="in"><code>sqlite3 survey.db &lt; gen-survey-database.sql</code></pre>


### Connecting and testing the SQLite DB installation


In order to connect to the created database, you need to start SQLite from within the folder where you have created the database. So, from within the "~/swc/sql" directory, type:


<pre class="in"><code>sqlite3 survey.db</code></pre>


The command "sqlite3 survey.db" opens the database itself and drops you into the database command line prompt. In SQLite a database is a flat file, which needs to be explicitly opened. SQLite is then started which is indicated by the change of the command line prompt to "sqlite", as shown in the following output:


<pre class="in"><code>/novice/sql$ sqlite3 survey.db 
SQLite version 3.7.15.2 2013-01-09 11:53:05
Enter &#34;.help&#34; for instructions
Enter SQL statements terminated with a &#34;;&#34;
sqlite&gt;  </code></pre>


Let us check the list the names and files of attached databases with the command ".databases", as shown in the following output:


<pre class="in"><code>sqlite&gt; .databases
seq  name             file                                                      
---  ---------------  ----------------------------------------------------------
0    main             ~/novice/sql/survey.db </code></pre>


You can check that the necessary tables "Person", "Survey", "Site" and "Visited" exist by typing:


<pre class="in"><code>.tables</code></pre>


and the output of ".tables" would look like this:


<pre class="in"><code>sqlite&gt; .tables
Person   Site     Survey   Visited</code></pre>


Now, you are done with the setup and you can proceed to the next lesson. You can conduct the following exercises in the current command line SQLite session. Since the IPython magic methods (the first rows in each command starting with %) will only work in an IPython notebook, you can omit them while using SQLite3 in the terminal. If you prefer to use the ipython notebook, you can quit SQLite3.


### How to exit the SQLite3 DB command line interface (CLI)


 To exit SQLite3, type:


<pre class="in"><code>sqlite&gt; .quit</code></pre>


### How to use the IPython notebook instead of the SQLite3 CLI


If you prefer to use an IPython notebook to follow the examples, check if IPython is installed on your machine. If it is not installed, follow the [installation instructions](http://ipython.org/install.html). 

If IPython is already installed on your local machine, to open a notebook, type "ipython notebook" from within the ~/swc/sql folder. 


<pre class="in"><code>~/swc/sql$ ipython notebook</code></pre>


The above command will launch the IPython kernel displaying an interactive notebook on your default browser, which can be edited as you learn. The commands shown in the next lessons can be executed as they are in the IPython notebook. When you are done, to keep your changes, remember to save your notebook.
