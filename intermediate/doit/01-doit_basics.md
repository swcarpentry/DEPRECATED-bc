---
layout: lesson
root: ../..
---

## Automating Tasks With "doit"


We're starting a project where we need to do some analysis of climate data. This analysis is going to require a number of steps, which all have to be carried out in the correct order. Our data is also updating all the time with new readings, so we don't want to have to keep track of which steps we have or have not remembered to re-run every time we update the source data.

In order to do this, we're going to use a python library called [doit](http://pydoit.org/). This lesson will cover the basics of doit, but doit has excellent [documentation](http://pydoit.org/contents.html) for those who are interested in more advanced usage.

### Objectives:


- Explain the difference between a dependency and a target
- Determine the order in which doit will execute a series of tasks
- Explain how automatic variables can reduce repetition in doit task definitions
- Write a simple doit task configuration file

## Basic doit files


Let's start by looking at the raw data we have to work with. There are two files containing data on monthly mean temperature and monthly total sunshine hours:


<pre class="in"><code>!ls *.txt</code></pre>

<div class="out"><pre class='out'><code>UK_Sunshine_data.txt  UK_Tmean_data.txt
</code></pre></div>


Now let's peek inside the mean temperatures file using head:


<pre class="in"><code>!head UK_Tmean_data.txt</code></pre>

<div class="out"><pre class='out'><code>UK Mean Temperature (Degrees C)
Areal series, starting from 1910
Allowances have been made for topographic, coastal and urban effects where relationships are found to exist.
Seasons: Winter=Dec-Feb, Spring=Mar-May, Summer=June-Aug, Autumn=Sept-Nov. (Winter: Year refers to Jan/Feb).
Monthly values are ranked and displayed to 1 dp and seasonal/annual values to 2 dp. Where values are equal, rankings are based in order of year descending.
Data are provisional from January 2012 &amp; Winter 2012 Last updated 01/12/2012

	JAN	Year	FEB	Year	MAR	Year	APR	Year	MAY	Year	JUN	Year	JUL	Year	AUG	Year	SEP	Year	OCT	Year	NOV	Year	DEC	Year	WIN	Year	SPR	Year	SUM	Year	AUT	Year	ANN	Year
	6.3	1916	6.8	1998	8	1938	10.7	2011	12.2	2008	15	1976	17.8	2006	17.3	1995	15.2	2006	12.2	2001	8.8	1994	6.9	1934	5.81	1989	9.15	2011	15.78	2006	11.39	2006	9.73	2006
	5.9	2007	5.9	1990	7.9	1957	10.2	2007	12	1992	14.9	1940	17.3	1983	17.1	1997	14.6	1949	11.8	1969	8.7	2011	6.6	1988	5.56	2007	9.05	2007	15.77	2003	11.26	2011	9.64	2011
</code></pre></div>


The data in this file is organized in a pretty terrible way. There are seven lines at the beginning of the file which explain the structure of the data, which is helpful (although it would be even better if they started with a comment character like #).

Essentially, there are two columns per month. The first column contains the mean temperature for that month and the second contains the year of the measurements. Every column is ordered by increasing temperature.

Thankfully, an old grad student left us a python script that can massage this data into a more useful format. Let's have a look at what this does:


<pre class="in"><code>!python reformat_weather_data.py UK_Tmean_data.txt | head</code></pre>

<div class="out"><pre class='out'><code>month,value
1910-01-01,2.6
1910-02-01,2.6
1910-03-01,4.0
1910-04-01,6.4
1910-05-01,9.5
1910-06-01,12.3
1910-07-01,14.0
1910-08-01,13.8
1910-09-01,11.8
</code></pre></div>


Much better. The first task in our analysis pipeline is to run this python script on the file `UK_Tmean_data.txt` and save it as a new file, `UK_Tmean_data.reformatted.txt`. We want to use the "doit" library for python to automatically perform this reformatting every time the raw data is updated.

First make sure doit is installed:


<pre class="in"><code>!pip install doit</code></pre>

<div class="out"><pre class='out'><code>Requirement already satisfied (use --upgrade to upgrade): doit in /usr/local/lib/python2.7/dist-packages
Requirement already satisfied (use --upgrade to upgrade): six in /usr/local/lib/python2.7/dist-packages (from doit)
Requirement already satisfied (use --upgrade to upgrade): pyinotify in /usr/lib/python2.7/dist-packages (from doit)
Cleaning up...
</code></pre></div>


Normally, we make a file containing the details of all our tasks inside of a python file. If you have had some experience of  make, this file is the equivalent of a makefile. If we then run the command `doit` in our terminal, doit will look for a configuration file called dodo.py in the current directory, read the tasks from the file and execute out those which are out of date. We can also use `doit -f <name_of_task_file.py>` to get doit to read a file which is not called dodo.py.

For the sake of convenience in this lesson, I'll be using some iPython magic to run doit code from the iPython notebook. In each case, the contents of the cell corresponds to what you would put in your `dodo.py` file.


<pre class="in"><code>%load_ext doitmagic</code></pre>


Here is our first doit file, containing just one task:


<pre class="in"><code>%%doit

# one_task.py

def task_reformat_temperature_data():
    &#34;&#34;&#34;Reformats the raw temperature data file for easier analysis&#34;&#34;&#34;
    
    return {
        &#39;file_dep&#39;: [&#39;UK_Tmean_data.txt&#39;],
        &#39;targets&#39;: [&#39;UK_Tmean_data.reformatted.txt&#39;],
        &#39;actions&#39;: [&#39;python reformat_weather_data.py UK_Tmean_data.txt &gt; UK_Tmean_data.reformatted.txt&#39;],
    }
</code></pre>

<div class="out"><pre class='out'><code>.  reformat_temperature_data
</code></pre></div>


The python function defines a single task that we want doit to carry out. All the function does is return a dictionary containing the configuration for this task. Lets look in more detail at the components of this configuration:

The task has one file dependency, or `file_dep` - this tells doit that the task depends on the `UK_Tmean_data.txt` file, so if that file has changed we need to re-run the task. 

It also has one `target` - this tells doit that the task creates the `UK_Tmean_data.reformatted.txt` file. If the `UK_Tmean_data.reformatted.txt` file doesn't exist, we need to run this task to create it. 

Finally, the task has one `action`. The `actions` part of the task definition is a list of commands to run when doit determines that the task is not up to date.

Now let's look at doit's output. Doit shows the name of each task on a seperate line, and since we only gave it one task we only get one line. Since we didn't explicitly give our task a name, doit guesses the name from the function name. The dot (`.`) before the task name means that doit determined that the task was actually run. We can run that cell again and see what changes.

Now the task name is preceded by two dashes (`--`), which means that doit found our task, but since the `UK_Tmean_data.reformatted.txt` file already exists and the `UK_Tmean_data.txt` file hasn't changed, it didn't run the task again.

We should check the new file to make sure that the task we wrote does what we want:


<pre class="in"><code>!head UK_Tmean_data.reformatted.txt</code></pre>

<div class="out"><pre class='out'><code>month,value
1910-01-01,2.6
1910-02-01,2.6
1910-03-01,4.0
1910-04-01,6.4
1910-05-01,9.5
1910-06-01,12.3
1910-07-01,14.0
1910-08-01,13.8
1910-09-01,11.8
</code></pre></div>


If we were only allowed one rule per file, this wouldn't be any simpler than typing commands by hand or putting them in little shell scripts. Luckily, doit allows us to put any number of rules in a single configuration file. 

Here is another doit file called two_tasks.py with rules to reformat both `UK_Tmean_data.txt` and `UK_Sunshine_data.txt`. These rules are identical except for the "Tmean" or "Sunshine" in the filenames; we'll see later how to combine these rules into one.


<pre class="in"><code>%%doit

# two_tasks.py

def task_reformat_temperature_data():
    &#34;&#34;&#34;Reformats the raw temperature data file for easier analysis&#34;&#34;&#34;
        
    return {
        &#39;file_dep&#39;: [&#39;UK_Tmean_data.txt&#39;],
        &#39;targets&#39;: [&#39;UK_Tmean_data.reformatted.txt&#39;],
        &#39;actions&#39;: [&#39;python reformat_weather_data.py UK_Tmean_data.txt &gt; UK_Tmean_data.reformatted.txt&#39;],
}

def task_reformat_sunshine_data():
    &#34;&#34;&#34;Reformats the raw sunshine data file for easier analysis&#34;&#34;&#34;

    return {
        &#39;file_dep&#39;: [&#39;UK_Sunshine_data.txt&#39;],
        &#39;targets&#39;: [&#39;UK_Sunshine_data.reformatted.txt&#39;],
        &#39;actions&#39;: [&#39;python reformat_weather_data.py UK_Sunshine_data.txt &gt; UK_Sunshine_data.reformatted.txt&#39;],
    }</code></pre>

<div class="out"><pre class='out'><code>-- reformat_temperature_data
.  reformat_sunshine_data
</code></pre></div>


Now we see that doit found both of our tasks. It determined that it didn't need to run the task that reformats the temperature data, but it did run our new task that reformats the sunshine data.

If we run the cell again, we should see that now doit decides it doesn't need to run either task.

One thing to note is that if there is no dependency to satisfy between the tasks then doit executes them in the order they are defined. It could also execute them in parallel if it had more than one processor to use - we'll return to this idea later.

Something else this example shows us is that a single thing can be a target in one rule, and a prerequisite in others. The dependencies between the files mentioned in the dodo.py make up a directed graph. In order for doit to run, this graph must not contain any cycles. For example, if X depends on Y, Y depends on Z, and Z depends on X, everything depends on something else, so there is nothing doit can execute first. If it detects a cycle in between tasks, doit will print an error message and stop.

As we noted earlier, there is a lot of redundancy in this file. Firstly, the file names are repeated in the task definition and the task's action. Luckily, doit gives us access to some variables when we are writing our tasks actions.

Doit uses python's `%` formatter to substitute a task's dependencies and targets in the string which defines the action. It works like this:


<pre class="in"><code>%%doit

# automatic_variables.py

def task_reformat_temperature_data():
    &#34;&#34;&#34;Reformats the raw temperature data file for easier analysis&#34;&#34;&#34;
    
    return {
        &#39;actions&#39;: [&#39;python reformat_weather_data.py %(dependencies)s &gt; %(targets)s&#39;],
        &#39;file_dep&#39;: [&#39;UK_Tmean_data.txt&#39;],
        &#39;targets&#39;: [&#39;UK_Tmean_data.reformatted.txt&#39;],
    }

def task_reformat_sunshine_data():
    &#34;&#34;&#34;Reformats the raw sunshine data file for easier analysis&#34;&#34;&#34;
    
    return {
        &#39;actions&#39;: [&#39;python reformat_weather_data.py %(dependencies)s &gt; %(targets)s&#39;],
        &#39;file_dep&#39;: [&#39;UK_Sunshine_data.txt&#39;],
        &#39;targets&#39;: [&#39;UK_Sunshine_data.reformatted.txt&#39;],
    }</code></pre>

<div class="out"><pre class='out'><code>-- reformat_temperature_data
-- reformat_sunshine_data
</code></pre></div>


This is better, but now the action is identical between the two tasks. Only the dependency and the target are different.

We'll remove the rest of the redundancy in the next section.

### Challenges:


1. Write a task that uses the unix "echo" command to create a new file called hello.txt, containing the text "Hello world!"
2. Given the following task configuration file, in what order would doit execute the tasks:


<pre class="in"><code>def task_giraffe():
    
    return {
            &#39;targets&#39; : [&#39;giraffe.txt&#39;],
            &#39;actions&#39; : [&#39;touch %(targets)s&#39;]
           }

def task_zebra():
    
    return {
            &#39;targets&#39; : [&#39;zebra.txt&#39;],
            &#39;file_dep&#39;: [&#39;lion.txt&#39;],
            &#39;actions&#39; : [&#39;touch %(targets)s&#39;]
           }

def task_lion():
    
    return {
            &#39;targets&#39; : [&#39;lion.txt&#39;],
            &#39;actions&#39; : [&#39;touch %(targets)s&#39;]
           }</code></pre>
