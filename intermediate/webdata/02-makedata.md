---
layout: lesson
root: ../..
---

## Creating Data


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Write Python programs that share data in a findable way.</li>
</ul>
</div>

### Step 1: Creating a Data Set


<div>
<p>In our previous lesson, we built functions called <code>get_annual_mean_temp_by_country</code> and <code>diff_records</code> to download temperature data for different countries and find annual differences. The next step is to share our findings with the world by making publishing the data sets we generate. To do this, we have to answer three questions:</p>
<ol style="list-style-type: decimal">
<li>How are we going to store the data?</li>
<li>How are people going to download it?</li>
<li>How are people going to <em>find</em> it?</li>
</ol>
<p>The first question is the easiest to answer: <code>diff_records</code> returns a list of (year, difference) pairs that we can write out as a CSV file:</p>
</div>


<div class="in">
<pre>import csv

def save_records(filename, records):
    &#39;&#39;&#39;Save a list of [year, temp] pairs as CSV.&#39;&#39;&#39;
    with open(filename, &#39;w&#39;) as raw:
        writer = csv.writer(raw)
        writer.writerows(records)</pre>
</div>


<div>
<p>Let's test it:</p>
</div>


<div class="in">
<pre>save_records(&#39;temp.csv&#39;, [[1, 2], [3, 4]])</pre>
</div>


<div class="in">
<pre>!cat temp.csv</pre>
</div>

<div class="out">
<pre>1,2
3,4
</pre>
</div>


<div class="challenges">
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Modify <code>save_records</code> so that it can be tested using <code>cStringIO</code>.</p></li>
<li><p>Should <code>save_records</code> check that every record in its input is the same length? Why or why not?</p></li>
</ol>
</div>

### Step 2: Publishing Data


<div>
<p>Now, where should this file go? The answer is clearly &quot;a server&quot;, since data on our laptop is only accessible when we're online (and probably not even then, since most people don't run a web server on their laptop). But where on the server, and what should we call it?</p>
<p>The answer to those questions depends on how the server is set up. On many multi-user Linux machines, users can create a directory called something like <code>public_html</code> under their home directory, and the web server will search in those directories. For example, if Nelle has a file called <code>thesis.pdf</code> in her <code>public_html</code> directory, the web server will find it when it gets the URL <code>http://the.server.name/u/nelle/thesis.pdf</code>. The specifics differ from one machine to the next, but the mechanism stays the same.</p>
<p>As for what we should call it, here we return to the key idea in REST: every data set should be identified by a &quot;guessable&quot; URL. In our case, we'll use the name <code>left-right.csv</code>, where <code>left</code> and <code>right</code> are the three-letter codes of the countries whose mean annual temperatures we are differencing. We can then tell people that if they want to compare Australia and Brazil, they should look for <code>http://the.server.name/u/nelle/AUS-BRA.csv</code>. (We use upper case to be consistent with the World Bank's API.)</p>
<p>But what's to prevent someone from creating a badly-named (and therefore unfindable) file? Someone could, for example, call <code>save_records('aus+bra.csv', records)</code>. To prevent this (or at least reduce the risk), let's modify <code>save_records</code> as follows:</p>
</div>


<div class="in">
<pre>import csv

def save_records(left, right, records):
    &#39;&#39;&#39;Save a list of [year, temp] pairs as CSV.&#39;&#39;&#39;
    filename = left + &#39;-&#39; + right + &#39;.csv&#39;
    with open(filename, &#39;w&#39;) as raw:
        writer = csv.writer(raw)
        writer.writerows(records)</pre>
</div>


<div>
<p>We can now call it like this:</p>
</div>


<div class="in">
<pre>save_records(&#39;AUS&#39;, &#39;BRA&#39;, [[1, 2], [3, 4]])</pre>
</div>


<div>
<p>and then check that the right output file has been created:</p>
</div>


<div class="in">
<pre>!cat AUS-BRA.csv</pre>
</div>

<div class="out">
<pre>1,2
3,4
</pre>
</div>


<div>
<p>Since we are bound to have the country codes anyway (having used them to look up our data), this is as little extra work as possible.</p>
</div>


<div class="challenges">
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li>Find out how to publish a file on your department's server.</li>
</ol>
</div>

### Step 3: Making Data Findable


<div>
<p>The final step in this lesson is to make the data we generate findable. It's not enough to tell people what the rule is for creating filenames, since that doesn't tell them what data sets we've actually generated. Instead, we need to create an <a href="../../gloss.html#index">index</a> to tell them what files exist. For reasons we will see in a moment, that index should also tell them when each data set was generated.</p>
<p>Here's the format we will use:</p>
</div>


<div class="in">
<pre>!cat index.csv</pre>
</div>

<div class="out">
<pre>2014-05-26,FRA,TCD,FRA-TCD.csv
2014-05-27,AUS,BRA,AUS-BRA.csv
2014-05-27,AUS,CAN,AUS-CAN.csv
2014-05-28,BRA,CAN,BRA-CAN.csv
</pre>
</div>


<div>
<p>The four columns in this file are self-explanatory, but why do we bother to include the filename? After all, we can re-generate it easily given the two country codes. The answer is that while <em>we</em> know the rule for generating filenames, other people's programs shouldn't have to. Explicit is better than implicit.</p>
<p>Here's a function that updates the index file every time we generate a new data file:</p>
</div>


<div class="in">
<pre>import time

def update_index(index_filename, left, right):
    &#39;&#39;&#39;Append a record to the index.&#39;&#39;&#39;

    # Read existing data.
    with open(index_filename, &#39;r&#39;) as raw:
        reader = csv.reader(raw)
        records = []
        for r in reader:
            records.append(r)
    
    # Create new record.
    timestamp = time.strftime(&#39;%Y-%m-%d&#39;)
    data_filename = left + &#39;-&#39; + right + &#39;.csv&#39;
    new_record = (timestamp, left, right, data_filename)
    
    # Save.
    records.append(new_record)
    with open(index_filename, &#39;w&#39;) as raw:
        writer = csv.writer(raw)
        writer.writerows(records)</pre>
</div>


<div>
<p>Let's test it:</p>
</div>


<div class="in">
<pre>!cat index.csv</pre>
</div>

<div class="out">
<pre>2014-05-26,FRA,TCD,FRA-TCD.csv
2014-05-27,AUS,BRA,AUS-BRA.csv
2014-05-27,AUS,CAN,AUS-CAN.csv
2014-05-28,BRA,CAN,BRA-CAN.csv
</pre>
</div>


<div class="in">
<pre>update_index(&#39;index.csv&#39;, &#39;TCD&#39;, &#39;CAN&#39;)</pre>
</div>


<div class="in">
<pre>!cat index.csv</pre>
</div>

<div class="out">
<pre>2014-05-26,FRA,TCD,FRA-TCD.csv
2014-05-27,AUS,BRA,AUS-BRA.csv
2014-05-27,AUS,CAN,AUS-CAN.csv
2014-05-28,BRA,CAN,BRA-CAN.csv
2014-05-29,TCD,CAN,TCD-CAN.csv
</pre>
</div>


<div class="challenges">
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Should <code>update_index</code> be called inside <code>save_records</code> so that the index is automatically updated every time a new data set is generated? Why or why not?</p></li>
<li><p><code>update_index</code> and <code>save_records</code> both construct the name of the data file. Refactor them to remove this redundancy.</p></li>
</ol>
</div>

### The Payoff


<div>
<p>Now that all of this is in place, it's easy for us—and other people—to do new and exciting things with our data. For example, we can easily write a small program that tells us what data sets are available that include information about a particular country <em>and</em> have been published since we last checked:</p>
</div>


<div class="in">
<pre>def what_is_available(index_file, country, after):
    &#39;&#39;&#39;What data files include a country and have been published since &#39;after&#39;?&#39;&#39;&#39;
    with open(index_file, &#39;r&#39;) as raw:
        reader = csv.reader(raw)
        filenames = []
        for record in reader:
            if (record[0] &gt;= after) and (record[1] == country or record[2] == country):
                filenames.append(record[3])
    return filenames

print what_is_available(&#39;index.csv&#39;, &#39;BRA&#39;, &#39;2014-05-27&#39;)</pre>
</div>

<div class="out">
<pre>[&#39;AUS-BRA.csv&#39;, &#39;BRA-CAN.csv&#39;]
</pre>
</div>


<div>
<p>This may not seem like a breakthrough, but it is actually an example of how the World Wide Web helps researchers do new kinds of science. With a little bit more work, we could create a file on <em>our</em> machine to record when we last ran <code>what_is_available</code> for each of several different sites that are producing data. Each time we run it, our program would:</p>
<ol style="list-style-type: decimal">
<li>read our local &quot;what to check&quot; file;</li>
<li>ask each data source whether it had any new data for us;</li>
<li>download and process that data; and</li>
<li>present us with a summary of the results.</li>
</ol>
<p>This is exactly how blogs work. Every blog reader keeps a list of blog URLs that it's supposed to check. When it is run, it goes to each of those sites and asks them for their index file (which is typically called something like <code>feed.xml</code>). It then checks the articles listed in that index against its local record of what has already been seen, then downloads any articles that are new.</p>
<p>By automating this process, blogging tools help us focus attention on things that are actually worth looking at.</p>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>Publish data by putting files with predictable names in a publicly-accessible location.</li>
<li>Create a machine-readable index to explicitly tell people what data sets are available.</li>
</ul>
</div>
