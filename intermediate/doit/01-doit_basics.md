---
layout: lesson
root: ../..
---

## Automating Tasks With "doit"


<div>
<p>We're starting a project where we need to do some analysis of climate data. This analysis is going to require a number of steps, which all have to be carried out in the correct order. Our data is also updating all the time with new readings, so we don't want to have to keep track of which steps we have or have not remembered to re-run every time we update the source data.</p>
<p>In order to do this, we're going to use a python library called <a href="http://pydoit.org/">doit</a>. This lesson will cover the basics of doit, but doit has excellent <a href="http://pydoit.org/contents.html">documentation</a> for those who are interested in more advanced usage.</p>
</div>

### Objectives:


<div>
<ul>
<li>Explain the difference between a dependency and a target</li>
<li>Determine the order in which doit will execute a series of tasks</li>
<li>Explain how automatic variables can reduce repetition in doit task definitions</li>
<li>Write a simple doit task configuration file</li>
</ul>
</div>

## Basic doit files


<div>
<p>Let's start by looking at the raw data we have to work with. There are two files containing data on monthly mean temperature and monthly total sunshine hours:</p>
</div>


<div class="in">
<pre>!ls *.txt</pre>
</div>

<div class="out">
<pre>UK_Sunshine_data.txt  UK_Tmean_data.txt
</pre>
</div>


<div>
<p>Now let's peek inside the mean temperatures file using head:</p>
</div>


<div class="in">
<pre>!head UK_Tmean_data.txt</pre>
</div>

<div class="out">
<pre>UK Mean Temperature (Degrees C)
Areal series, starting from 1910
Allowances have been made for topographic, coastal and urban effects where relationships are found to exist.
Seasons: Winter=Dec-Feb, Spring=Mar-May, Summer=June-Aug, Autumn=Sept-Nov. (Winter: Year refers to Jan/Feb).
Monthly values are ranked and displayed to 1 dp and seasonal/annual values to 2 dp. Where values are equal, rankings are based in order of year descending.
Data are provisional from January 2012 &amp; Winter 2012 Last updated 01/12/2012

	JAN	Year	FEB	Year	MAR	Year	APR	Year	MAY	Year	JUN	Year	JUL	Year	AUG	Year	SEP	Year	OCT	Year	NOV	Year	DEC	Year	WIN	Year	SPR	Year	SUM	Year	AUT	Year	ANN	Year
	6.3	1916	6.8	1998	8	1938	10.7	2011	12.2	2008	15	1976	17.8	2006	17.3	1995	15.2	2006	12.2	2001	8.8	1994	6.9	1934	5.81	1989	9.15	2011	15.78	2006	11.39	2006	9.73	2006
	5.9	2007	5.9	1990	7.9	1957	10.2	2007	12	1992	14.9	1940	17.3	1983	17.1	1997	14.6	1949	11.8	1969	8.7	2011	6.6	1988	5.56	2007	9.05	2007	15.77	2003	11.26	2011	9.64	2011
</pre>
</div>


<div>
<p>The data in this file is organized in a pretty terrible way. There are seven lines at the beginning of the file which explain the structure of the data, which is helpful (although it would be even better if they started with a comment character like #).</p>
<p>Essentially, there are two columns per month. The first column contains the mean temperature for that month and the second contains the year of the measurements. Every column is ordered by increasing temperature.</p>
<p>Thankfully, an old grad student left us a python script that can massage this data into a more useful format. Let's have a look at what this does:</p>
</div>


<div class="in">
<pre>!python reformat_weather_data.py UK_Tmean_data.txt | head</pre>
</div>

<div class="out">
<pre>month,value
1910-01-01,2.6
1910-02-01,2.6
1910-03-01,4.0
1910-04-01,6.4
1910-05-01,9.5
1910-06-01,12.3
1910-07-01,14.0
1910-08-01,13.8
1910-09-01,11.8
</pre>
</div>


<div>
<p>Much better. The first task in our analysis pipeline is to run this python script on the file <code>UK_Tmean_data.txt</code> and save it as a new file, <code>UK_Tmean_data.reformatted.txt</code>. We want to use the &quot;doit&quot; library for python to automatically perform this reformatting every time the raw data is updated.</p>
<p>First make sure doit is installed:</p>
</div>


<div class="in">
<pre>!pip install doit</pre>
</div>

<div class="out">
<pre>Requirement already satisfied (use --upgrade to upgrade): doit in /usr/local/lib/python2.7/dist-packages
Requirement already satisfied (use --upgrade to upgrade): six in /usr/local/lib/python2.7/dist-packages (from doit)
Requirement already satisfied (use --upgrade to upgrade): pyinotify in /usr/lib/python2.7/dist-packages (from doit)
Cleaning up...
</pre>
</div>


<div>
<p>Normally, we make a file containing the details of all our tasks inside of a python file. If you have had some experience of make, this file is the equivalent of a makefile. If we then run the command <code>doit</code> in our terminal, doit will look for a configuration file called dodo.py in the current directory, read the tasks from the file and execute out those which are out of date. We can also use <code>doit -f &lt;name_of_task_file.py&gt;</code> to get doit to read a file which is not called dodo.py.</p>
<p>For the sake of convenience in this lesson, I'll be using some iPython magic to run doit code from the iPython notebook. In each case, the contents of the cell corresponds to what you would put in your <code>dodo.py</code> file.</p>
</div>


<div class="in">
<pre>%load_ext doitmagic</pre>
</div>


<div>
<p>Here is our first doit file, containing just one task:</p>
</div>


<div class="in">
<pre>%%doit

# one_task.py

def task_reformat_temperature_data():
    &#34;&#34;&#34;Reformats the raw temperature data file for easier analysis&#34;&#34;&#34;
    
    return {
        &#39;file_dep&#39;: [&#39;UK_Tmean_data.txt&#39;],
        &#39;targets&#39;: [&#39;UK_Tmean_data.reformatted.txt&#39;],
        &#39;actions&#39;: [&#39;python reformat_weather_data.py UK_Tmean_data.txt &gt; UK_Tmean_data.reformatted.txt&#39;],
    }
</pre>
</div>

<div class="out">
<pre>.  reformat_temperature_data
</pre>
</div>


<div>
<p>The python function defines a single task that we want doit to carry out. All the function does is return a dictionary containing the configuration for this task. Lets look in more detail at the components of this configuration:</p>
<p>The task has one file dependency, or <code>file_dep</code> - this tells doit that the task depends on the <code>UK_Tmean_data.txt</code> file, so if that file has changed we need to re-run the task.</p>
<p>It also has one <code>target</code> - this tells doit that the task creates the <code>UK_Tmean_data.reformatted.txt</code> file. If the <code>UK_Tmean_data.reformatted.txt</code> file doesn't exist, we need to run this task to create it.</p>
<p>Finally, the task has one <code>action</code>. The <code>actions</code> part of the task definition is a list of commands to run when doit determines that the task is not up to date.</p>
<p>Now let's look at doit's output. Doit shows the name of each task on a seperate line, and since we only gave it one task we only get one line. Since we didn't explicitly give our task a name, doit guesses the name from the function name. The dot (<code>.</code>) before the task name means that doit determined that the task was actually run. We can run that cell again and see what changes.</p>
<p>Now the task name is preceded by two dashes (<code>--</code>), which means that doit found our task, but since the <code>UK_Tmean_data.reformatted.txt</code> file already exists and the <code>UK_Tmean_data.txt</code> file hasn't changed, it didn't run the task again.</p>
<p>We should check the new file to make sure that the task we wrote does what we want:</p>
</div>


<div class="in">
<pre>!head UK_Tmean_data.reformatted.txt</pre>
</div>

<div class="out">
<pre>month,value
1910-01-01,2.6
1910-02-01,2.6
1910-03-01,4.0
1910-04-01,6.4
1910-05-01,9.5
1910-06-01,12.3
1910-07-01,14.0
1910-08-01,13.8
1910-09-01,11.8
</pre>
</div>


<div>
<p>If we were only allowed one rule per file, this wouldn't be any simpler than typing commands by hand or putting them in little shell scripts. Luckily, doit allows us to put any number of rules in a single configuration file.</p>
<p>Here is another doit file called two_tasks.py with rules to reformat both <code>UK_Tmean_data.txt</code> and <code>UK_Sunshine_data.txt</code>. These rules are identical except for the &quot;Tmean&quot; or &quot;Sunshine&quot; in the filenames; we'll see later how to combine these rules into one.</p>
</div>


<div class="in">
<pre>%%doit

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
    }</pre>
</div>

<div class="out">
<pre>-- reformat_temperature_data
.  reformat_sunshine_data
</pre>
</div>


<div>
<p>Now we see that doit found both of our tasks. It determined that it didn't need to run the task that reformats the temperature data, but it did run our new task that reformats the sunshine data.</p>
<p>If we run the cell again, we should see that now doit decides it doesn't need to run either task.</p>
<p>One thing to note is that if there is no dependency to satisfy between the tasks then doit executes them in the order they are defined. It could also execute them in parallel if it had more than one processor to use - we'll return to this idea later.</p>
<p>Something else this example shows us is that a single thing can be a target in one rule, and a prerequisite in others. The dependencies between the files mentioned in the dodo.py make up a directed graph. In order for doit to run, this graph must not contain any cycles. For example, if X depends on Y, Y depends on Z, and Z depends on X, everything depends on something else, so there is nothing doit can execute first. If it detects a cycle in between tasks, doit will print an error message and stop.</p>
<p>As we noted earlier, there is a lot of redundancy in this file. Firstly, the file names are repeated in the task definition and the task's action. Luckily, doit gives us access to some variables when we are writing our tasks actions.</p>
<p>Doit uses python's <code>%</code> formatter to substitute a task's dependencies and targets in the string which defines the action. It works like this:</p>
</div>


<div class="in">
<pre>%%doit

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
    }</pre>
</div>

<div class="out">
<pre>-- reformat_temperature_data
-- reformat_sunshine_data
</pre>
</div>


<div>
<p>This is better, but now the action is identical between the two tasks. Only the dependency and the target are different.</p>
<p>We'll remove the rest of the redundancy in the next section.</p>
</div>

### Challenges:


<div>
<ol style="list-style-type: decimal">
<li>Write a task that uses the unix &quot;echo&quot; command to create a new file called hello.txt, containing the text &quot;Hello world!&quot;</li>
<li>Given the following task configuration file, in what order would doit execute the tasks:</li>
</ol>
</div>


<div class="in">
<pre>def task_giraffe():
    
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
           }</pre>
</div>
