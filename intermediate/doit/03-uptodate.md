---
layout: lesson
root: ../..
---

## Determining if doit needs to run a task

### Objectives:


<div>
<ul>
<li>Explain how doit decides if a task is up to date when that task depends on another file</li>
<li>Explain how this is decided if a task does not have any dependencies</li>
<li>Explain how we configure a task to change the way doit will decide if it is up to date</li>
<li>Explain how doit decides which functions are tasks, and which are not</li>
</ul>
</div>


<div>
<p>Here is one version of the script, which now downloads both the raw data files. The download_data task was a bit long, so I refactored the part which calculates all the correct file names into a new python function.</p>
</div>


<div class="in">
<pre>%load_ext doitmagic</pre>
</div>


<div class="in">
<pre>%%doit

# download_all_data.py

data_sets = [&#39;Tmean&#39;, &#39;Sunshine&#39;]

def get_data_file_parameters(data_type):
    &#34;&#34;&#34;Takes a string describing the type of climate data, returns url and file name for that data&#34;&#34;&#34;
    
    base_url = &#39;http://www.metoffice.gov.uk/climate/uk/datasets/{0}/ranked/UK.txt&#39;
    data_url = base_url.format(data_type)
    data_target = &#39;UK_{0}_data.txt&#39;.format(data_type)
    return data_url, data_target

def task_download_data():
    &#34;&#34;&#34;Downloads all raw data files from the Met Office website&#34;&#34;&#34;

    for data_type in data_sets:
        data_url, data_target = get_data_file_parameters(data_type)
        yield {
            &#39;actions&#39;: [&#39;wget -O %(targets)s {0}&#39;.format(data_url)],
            &#39;targets&#39;: [ data_target ],
            &#39;name&#39; : data_type,
        }

def task_reformat_data():
    &#34;&#34;&#34;Reformats all raw files for easier analysis&#34;&#34;&#34;

    for data_type in data_sets:
        yield {
            &#39;actions&#39;: [&#39;python reformat_weather_data.py %(dependencies)s &gt; %(targets)s&#39;],
            &#39;file_dep&#39;: [&#39;UK_{}_data.txt&#39;.format(data_type)],
            &#39;targets&#39;: [&#39;UK_{}_data.reformatted.txt&#39;.format(data_type)],
            &#39;name&#39;: &#39;UK_{}_data.txt&#39;.format(data_type),
        }</pre>
</div>

<div class="out">
<pre>.  download_data:Tmean
.  download_data:Sunshine
-- reformat_data:UK_Sunshine_data.txt
-- reformat_data:UK_Tmean_data.txt
--2014-04-05 12:08:57--  http://www.metoffice.gov.uk/climate/uk/datasets/Tmean/ranked/UK.txt
Resolving www.metoffice.gov.uk (www.metoffice.gov.uk)... 23.63.99.234, 23.63.99.216
Connecting to www.metoffice.gov.uk (www.metoffice.gov.uk)|23.63.99.234|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 25576 (25K) [text/plain]
Saving to: ‘UK_Tmean_data.txt’

     0K .......... .......... ....                            100% 1.51M=0.02s

2014-04-05 12:08:57 (1.51 MB/s) - ‘UK_Tmean_data.txt’ saved [25576/25576]

--2014-04-05 12:08:57--  http://www.metoffice.gov.uk/climate/uk/datasets/Sunshine/ranked/UK.txt
Resolving www.metoffice.gov.uk (www.metoffice.gov.uk)... 23.63.99.216, 23.63.99.234
Connecting to www.metoffice.gov.uk (www.metoffice.gov.uk)|23.63.99.216|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 20986 (20K) [text/plain]
Saving to: ‘UK_Sunshine_data.txt’

     0K .......... ..........                                 100% 6.48M=0.003s

2014-04-05 12:08:57 (6.48 MB/s) - ‘UK_Sunshine_data.txt’ saved [20986/20986]

</pre>
</div>


<div>
<p>There are two things to notice here. Firstly, my new python function <code>get_data_file_parameters</code> doesn't start with <code>task_</code>, so doit doesn't try to run it as a task.</p>
<p>Secondly, no matter how many times we run this script, doit always re-downloads the data files. Since they don't depend on anything, doit doesn't know how to check that they are up to date. This is why it always re-downloads them.</p>
<p>If we were trying to do all these tasks ourselves, we would probably want to update our monthly temperature data every month. Doit lets us define another parameter in our task configuration dictionary, called <code>uptodate</code>. This should define a python function which will tell doit whether our task needs to be re-run.</p>
<p>This is one of the big advantages of doit: any python you can write is valid. In fact, even using task generators to create sub-tasks is just a convention. Any python script you can write to make your task configuration dictionaries can be made to work with doit.</p>
<p>Back to our problem about keeping our raw data files fresh. You could write your own function that checks how old the raw data files are, but thankfully doit comes with a utility function for doing just this. Lets set our data files to expire after four weeks:</p>
</div>


<div class="in">
<pre>%%doit

# monthly_raw_data_update.py

import datetime
from doit.tools import timeout 

data_sets = [&#39;Tmean&#39;, &#39;Sunshine&#39;]

def get_data_file_parameters(data_type):
    &#34;&#34;&#34;Takes a string describing the type of climate data, returns url and file name for that data&#34;&#34;&#34;

    base_url = &#39;http://www.metoffice.gov.uk/climate/uk/datasets/{0}/ranked/UK.txt&#39;
    data_url = base_url.format(data_type)
    data_target = &#39;UK_{0}_data.txt&#39;.format(data_type)
    return data_url, data_target

def task_download_data():
    &#34;&#34;&#34;Downloads all raw data files from the Met Office website&#34;&#34;&#34;

    for data_type in data_sets:
        data_url, data_target = get_data_file_parameters(data_type)
        yield {
            &#39;actions&#39;: [&#39;wget -O %(targets)s {0}&#39;.format(data_url)],
            &#39;targets&#39;: [ data_target ],
            &#39;name&#39; : data_type,
            &#39;uptodate&#39;: [timeout(datetime.timedelta(weeks=4))],

        }

def task_reformat_data():
    &#34;&#34;&#34;Reformats all raw files for easier analysis&#34;&#34;&#34;

    for data_type in data_sets:
        yield {
            &#39;actions&#39;: [&#39;python reformat_weather_data.py %(dependencies)s &gt; %(targets)s&#39;],
            &#39;file_dep&#39;: [&#39;UK_{}_data.txt&#39;.format(data_type)],
            &#39;targets&#39;: [&#39;UK_{}_data.reformatted.txt&#39;.format(data_type)],
            &#39;name&#39;: &#39;UK_{}_data.txt&#39;.format(data_type),
        }</pre>
</div>

<div class="out">
<pre>.  download_data:Tmean
.  download_data:Sunshine
-- reformat_data:UK_Sunshine_data.txt
-- reformat_data:UK_Tmean_data.txt
--2014-04-05 12:08:57--  http://www.metoffice.gov.uk/climate/uk/datasets/Tmean/ranked/UK.txt
Resolving www.metoffice.gov.uk (www.metoffice.gov.uk)... 23.63.99.234, 23.63.99.216
Connecting to www.metoffice.gov.uk (www.metoffice.gov.uk)|23.63.99.234|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 25576 (25K) [text/plain]
Saving to: ‘UK_Tmean_data.txt’

     0K .......... .......... ....                            100% 1.87M=0.01s

2014-04-05 12:08:57 (1.87 MB/s) - ‘UK_Tmean_data.txt’ saved [25576/25576]

--2014-04-05 12:08:57--  http://www.metoffice.gov.uk/climate/uk/datasets/Sunshine/ranked/UK.txt
Resolving www.metoffice.gov.uk (www.metoffice.gov.uk)... 23.63.99.216, 23.63.99.234
Connecting to www.metoffice.gov.uk (www.metoffice.gov.uk)|23.63.99.216|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 20986 (20K) [text/plain]
Saving to: ‘UK_Sunshine_data.txt’

     0K .......... ..........                                 100% 6.86M=0.003s

2014-04-05 12:08:57 (6.86 MB/s) - ‘UK_Sunshine_data.txt’ saved [20986/20986]

</pre>
</div>


<div>
<p>When we run this new script for the second time, doit knows that the raw data files are less than four weeks old, so it doesn't download them again.</p>
<p>The big advantage of defining all our tasks in this way, is that it becomes much easier to add a new dataset. Lets download some rainfall data by adding <code>Rainfall</code> to our list of datasets:</p>
</div>


<div class="in">
<pre>%%doit

# rainfall_data.py

import datetime
from doit.tools import timeout 

data_sets = [&#39;Tmean&#39;, &#39;Sunshine&#39;, &#39;Rainfall&#39;]

def get_data_file_parameters(data_type):
    &#34;&#34;&#34;Takes a string describing the type of climate data, returns url and file name for that data&#34;&#34;&#34;

    base_url = &#39;http://www.metoffice.gov.uk/climate/uk/datasets/{0}/ranked/UK.txt&#39;
    data_url = base_url.format(data_type)
    data_target = &#39;UK_{0}_data.txt&#39;.format(data_type)
    return data_url, data_target

def task_download_data():
    &#34;&#34;&#34;Downloads all raw data files from the Met Office website&#34;&#34;&#34;

    for data_type in data_sets:
        data_url, data_target = get_data_file_parameters(data_type)
        yield {
            &#39;actions&#39;: [&#39;wget -O %(targets)s {0}&#39;.format(data_url)],
            &#39;targets&#39;: [ data_target ],
            &#39;name&#39; : data_type,
            &#39;uptodate&#39;: [timeout(datetime.timedelta(weeks=4))],

        }

def task_reformat_data():
    &#34;&#34;&#34;Reformats all raw files for easier analysis&#34;&#34;&#34;

    for data_type in data_sets:
        yield {
            &#39;actions&#39;: [&#39;python reformat_weather_data.py %(dependencies)s &gt; %(targets)s&#39;],
            &#39;file_dep&#39;: [&#39;UK_{}_data.txt&#39;.format(data_type)],
            &#39;targets&#39;: [&#39;UK_{}_data.reformatted.txt&#39;.format(data_type)],
            &#39;name&#39;: &#39;UK_{}_data.txt&#39;.format(data_type),
        }</pre>
</div>

<div class="out">
<pre>.  download_data:Rainfall
-- download_data:Tmean
-- download_data:Sunshine
.  reformat_data:UK_Rainfall_data.txt
-- reformat_data:UK_Sunshine_data.txt
-- reformat_data:UK_Tmean_data.txt
--2014-04-05 12:08:58--  http://www.metoffice.gov.uk/climate/uk/datasets/Rainfall/ranked/UK.txt
Resolving www.metoffice.gov.uk (www.metoffice.gov.uk)... 23.63.99.234, 23.63.99.216
Connecting to www.metoffice.gov.uk (www.metoffice.gov.uk)|23.63.99.234|:80... connected.
HTTP request sent, awaiting response... 200 OK
Length: 25518 (25K) [text/plain]
Saving to: ‘UK_Rainfall_data.txt’

     0K .......... .......... ....                            100% 1.68M=0.01s

2014-04-05 12:08:58 (1.68 MB/s) - ‘UK_Rainfall_data.txt’ saved [25518/25518]

</pre>
</div>


<div>
<p>Now we have some shiny new rainfall data, properly formatted for further analysis:</p>
</div>


<div class="in">
<pre>!head UK_Rainfall_data.reformatted.txt</pre>
</div>

<div class="out">
<pre>month,value
1910-01-01,111.4
1910-02-01,79.5
1910-03-01,75.5
1910-04-01,69.1
1910-05-01,66.4
1910-06-01,65.3
1910-07-01,81.6
1910-08-01,90.3
1910-09-01,92.0
</pre>
</div>

### Challenges:


<div>
<ul>
<li>Edit the rainfall_data.py file so that doit downloads the raw data files if they are older than 20 seconds</li>
<li>What will happen if we move the 'uptodate' configuration line from the download_data task to the reformat_data task?</li>
<li>Try it out and see if you were right!</li>
</ul>
</div>
