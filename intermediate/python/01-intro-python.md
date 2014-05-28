---
layout: lesson
root: ../..
---

# Introduction


<div>
<p>This material assumes that you have programmed before. This first lecture provides a quick introduction to programming in Python for those who either haven't used Python before or need a quick refresher.</p>
<p>Let's start with a hypothetical problem we want to solve. We are interested in understanding the relationship between the weather and the number of mosquitos occuring in a particular year so that we can plan mosquito control measures accordingly. Since we want to apply these mosquito control measures at a number of different sites we need to understand both the relationship at a particular site and whether or not it is consistent across sites. The data we have to address this problem comes from the local government and are stored in tables in comma-separated values (CSV) files. Each file holds the data for a single location, each row holds the information for a single year at that location, and the columns hold the data on both mosquito numbers and the average temperature and rainfall from the beginning of mosquito breeding season. The first few rows of our first file look like:</p>
<pre><code>year,temperature,rainfall,mosquitos
2001,87,222,198
2002,72,103,105
2003,77,176,166</code></pre>
</div>

## Objectives


<div>
<ul>
<li>Conduct variable assignment, looping, and conditionals in Python</li>
<li>Use an external Python library</li>
<li>Read tabular data from a file</li>
<li>Subset and perform analysis on data</li>
<li>Display simple graphs</li>
</ul>
</div>

## Loading Data


<div>
<p>In order to load the data, we need to import a library called Pandas that knows how to operate on tables of data.</p>
</div>


<div class="in">
<pre>import pandas</pre>
</div>


<div>
<p>We can now use Pandas to read our data file.</p>
</div>


<div class="in">
<pre>pandas.read_csv(&#39;mosquito_data_A1.csv&#39;)</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
0  2001           87       222        198
1  2002           72       103        105
2  2003           77       176        166
3  2004           89       236        210
4  2005           88       283        242
5  2006           89       151        147
6  2007           71       121        117
7  2008           88       267        232
8  2009           85       211        191
9  2010           75       101        106</pre>
</div>


<div>
<p>The <code>read_csv()</code> function belongs to the <code>pandas</code> library. In order to run it we need to tell Python that it is part of <code>pandas</code> and we do this using the dot notation, which is used everywhere in Python to refer to parts of larger things.</p>
<p>When we are finished typing and press Shift+Enter, the notebook runs our command and shows us its output. In this case, the output is the data we just loaded.</p>
<p>Our call to <code>pandas.read_csv()</code> read data into memory, but didn't save it anywhere. To do that, we need to assign the array to a variable. In Python we use <code>=</code> to assign a new value to a variable like this:</p>
</div>


<div class="in">
<pre>data = pandas.read_csv(&#39;mosquito_data_A1.csv&#39;)</pre>
</div>


<div>
<p>This statement doesn't produce any output because assignment doesn't display anything. If we want to check that our data has been loaded, we can print the variable's value:</p>
</div>


<div class="in">
<pre>print data</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
0  2001           87       222        198
1  2002           72       103        105
2  2003           77       176        166
3  2004           89       236        210
4  2005           88       283        242
5  2006           89       151        147
6  2007           71       121        117
7  2008           88       267        232
8  2009           85       211        191
9  2010           75       101        106
</pre>
</div>


<div>
<p><code>print data</code> tells Python to display the text. Alternatively we could just include <code>data</code> as the last value in a code cell:</p>
</div>


<div class="in">
<pre>data</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
0  2001           87       222        198
1  2002           72       103        105
2  2003           77       176        166
3  2004           89       236        210
4  2005           88       283        242
5  2006           89       151        147
6  2007           71       121        117
7  2008           88       267        232
8  2009           85       211        191
9  2010           75       101        106</pre>
</div>


<div>
<p>This tells the IPython Notebook to display the <code>data</code> object, which is why we see a pretty formated table.</p>
</div>

## Manipulating data


<div>
<p>Once we have imported the data we can start doing things with it. First, let's ask what type of thing <code>data</code> refers to:</p>
</div>


<div class="in">
<pre>print type(data)</pre>
</div>

<div class="out">
<pre>&lt;class &#39;pandas.core.frame.DataFrame&#39;&gt;
</pre>
</div>


<div>
<p>The data is stored in a data structure called a DataFrame. There are other kinds of data structures that are also commonly used in scientific computing including Numpy arrays, and Numpy matrices, which can be used for doing linear algebra.</p>
</div>


<div>
<p>We can select an individual column of data using its name:</p>
</div>


<div class="in">
<pre>print data[&#39;year&#39;]</pre>
</div>

<div class="out">
<pre>0    2001
1    2002
2    2003
3    2004
4    2005
5    2006
6    2007
7    2008
8    2009
9    2010
Name: year, dtype: int64
</pre>
</div>


<div>
<p>Or we can select several columns of data at once:</p>
</div>


<div class="in">
<pre>print data[[&#39;rainfall&#39;, &#39;temperature&#39;]]</pre>
</div>

<div class="out">
<pre>   rainfall  temperature
0       222           87
1       103           72
2       176           77
3       236           89
4       283           88
5       151           89
6       121           71
7       267           88
8       211           85
9       101           75
</pre>
</div>


<div>
<p>We can also select subsets of rows using slicing. Say we just want the first two rows of data:</p>
</div>


<div class="in">
<pre>print data[0:2]</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
0  2001           87       222        198
1  2002           72       103        105
</pre>
</div>


<div>
<p>There are a couple of important things to note here. First, Python indexing starts at zero. In contrast, programming languages like R and MATLAB start counting at 1, because that's what human beings have done for thousands of years. Languages in the C family (including C++, Java, Perl, and Python) count from 0 because that's simpler for computers to do. This means that if we have 5 things in Python they are numbered 0, 1, 2, 3, 4, and the first row in a data frame is always row 0.</p>
<p>The other thing to note is that the subset of rows starts at the first value and goes up to, but does not include, the second value. Again, the up-to-but-not-including takes a bit of getting used to, but the rule is that the difference between the upper and lower bounds is the number of values in the slice.</p>
</div>


<div>
<p>One thing that we can't do with this syntax is directly ask for the data from a single row:</p>
</div>


<div class="in">
<pre>data[1]</pre>
</div>

<div class="out">
<pre>---------------------------------------------------------------------------
KeyError                                  Traceback (most recent call last)
&lt;ipython-input-10-c805864c0d75&gt; in &lt;module&gt;()
----&gt; 1 data[1]

/usr/lib/python2.7/dist-packages/pandas/core/frame.pyc in __getitem__(self, key)
   2001             # get column
   2002             if self.columns.is_unique:
-&gt; 2003                 return self._get_item_cache(key)
   2004 
   2005             # duplicate columns

/usr/lib/python2.7/dist-packages/pandas/core/generic.pyc in _get_item_cache(self, item)
    665             return cache[item]
    666         except Exception:
--&gt; 667             values = self._data.get(item)
    668             res = self._box_item_values(item, values)
    669             cache[item] = res

/usr/lib/python2.7/dist-packages/pandas/core/internals.pyc in get(self, item)
   1653     def get(self, item):
   1654         if self.items.is_unique:
-&gt; 1655             _, block = self._find_block(item)
   1656             return block.get(item)
   1657         else:

/usr/lib/python2.7/dist-packages/pandas/core/internals.pyc in _find_block(self, item)
   1933 
   1934     def _find_block(self, item):
-&gt; 1935         self._check_have(item)
   1936         for i, block in enumerate(self.blocks):
   1937             if item in block:

/usr/lib/python2.7/dist-packages/pandas/core/internals.pyc in _check_have(self, item)
   1940     def _check_have(self, item):
   1941         if item not in self.items:
-&gt; 1942             raise KeyError(&#39;no item named %s&#39; % com.pprint_thing(item))
   1943 
   1944     def reindex_axis(self, new_axis, method=None, axis=0, copy=True):

KeyError: u&#39;no item named 1&#39;</pre>
</div>


<div>
<p>This is because there are several things that we could mean by <code>data[1]</code> so if we want a single row we can either take a slice that returns a single row:</p>
</div>


<div class="in">
<pre>print data[1:2]</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
1  2002           72       103        105
</pre>
</div>


<div>
<p>or use the <code>.iloc</code> method, which stands for &quot;integer location&quot; since we are looking up the row based on its integer index.</p>
</div>


<div class="in">
<pre>print data.iloc[1]</pre>
</div>

<div class="out">
<pre>year           2002
temperature      72
rainfall        103
mosquitos       105
Name: 1, dtype: int64
</pre>
</div>


<div>
<p>We can also use this same syntax for getting larger subsets of rows:</p>
</div>


<div class="in">
<pre>print data.iloc[1:3]</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
1  2002           72       103        105
2  2003           77       176        166
</pre>
</div>


<div>
<p>We can also subset the data based on the value of other rows:</p>
</div>


<div class="in">
<pre>print data[&#39;temperature&#39;][data[&#39;year&#39;] &gt; 2005]</pre>
</div>

<div class="out">
<pre>5    89
6    71
7    88
8    85
9    75
Name: temperature, dtype: int64
</pre>
</div>


<div>
<p>Data frames also know how to perform common mathematical operations on their values. If we want to find the average value for each variable, we can just ask the data frame for its mean values</p>
</div>


<div class="in">
<pre>print data.mean()</pre>
</div>

<div class="out">
<pre>year           2005.5
temperature      82.1
rainfall        187.1
mosquitos       171.4
dtype: float64
</pre>
</div>


<div>
<p>Data frames have lots of useful methods:</p>
</div>


<div class="in">
<pre>print data.max()</pre>
</div>

<div class="out">
<pre>year           2010
temperature      89
rainfall        283
mosquitos       242
dtype: int64
</pre>
</div>


<div class="in">
<pre>print data[&#39;temperature&#39;].min()</pre>
</div>

<div class="out">
<pre>71
</pre>
</div>


<div class="in">
<pre>print data[&#39;mosquitos&#39;][1:3].std()</pre>
</div>

<div class="out">
<pre>43.1335136524
</pre>
</div>

### Challenge


<div>
<p>Import the data from <code>mosquito_data_A2.csv</code>, create a new variable that holds a data frame with only the weather data, and print the means and standard deviations for the weather variables.</p>
</div>

## Loops


<div>
<p>Once we have some data we often want to be able to loop over it to perform the same operation repeatedly. A <code>for</code> loop in Python takes the general form</p>
<pre><code>for item in list:
    do_something</code></pre>
<p>So if we want to loop over the temperatures and print out there values in degrees Celcius (instead of Farenheit) we can use:</p>
</div>


<div class="in">
<pre>temps = data[&#39;temperature&#39;]
for temp_in_f in temps:
    temp_in_c = (temp_in_f - 32) * 5 / 9.0
    print temp_in_c</pre>
</div>

<div class="out">
<pre>30.5555555556
22.2222222222
25.0
31.6666666667
31.1111111111
31.6666666667
21.6666666667
31.1111111111
29.4444444444
23.8888888889
</pre>
</div>


<div>
<p>That looks good, but why did we use 9.0 instead of 9? The reason is that computers store integers and numbers with decimals as different types: integers and floating point numbers (or floats). Addition, subtraction and multiplication work on both as we'd expect, but division works differently. If we divide one integer by another, we get the quotient without the remainder:</p>
</div>


<div class="in">
<pre>print &#39;10/3 is:&#39;, 10 / 3</pre>
</div>

<div class="out">
<pre>10/3 is: 3
</pre>
</div>


<div>
<p>If either part of the division is a float, on the other hand, the computer creates a floating-point answer:</p>
</div>


<div class="in">
<pre>print &#39;10/3.0 is:&#39;, 10 / 3.0</pre>
</div>

<div class="out">
<pre>10/3.0 is: 3.33333333333
</pre>
</div>


<div>
<p>The computer does this for historical reasons: integer operations were much faster on early machines, and this behavior is actually useful in a lot of situations. However, it's still confusing, so Python 3 produces a floating-point answer when dividing integers if it needs to. We're still using Python 2.7 in this class, so if we want 5/9 to give us the right answer, we have to write it as 5.0/9, 5/9.0, or some other variation.</p>
</div>

## Conditionals


<div>
<p>The other standard thing we need to know how to do in Python is conditionals, or if/then/else statements. In Python the basic syntax is:</p>
<pre class="sourceCode python"><code class="sourceCode python"><span class="kw">if</span> condition:
    do_something</code></pre>
<p>So if we want to loop over the temperatures and print out only those temperatures that are greater than 80 degrees we would use:</p>
</div>


<div class="in">
<pre>temp = data[&#39;temperature&#39;][0]
if temp &gt; 80:
    print &#34;The temperature is greater than 80&#34;</pre>
</div>

<div class="out">
<pre>The temperature is greater than 80
</pre>
</div>


<div>
<p>We can also use <code>==</code> for equality, <code>&lt;=</code> for less than or equal to, <code>&gt;=</code> for greater than or equal to, and <code>!=</code> for not equal to.</p>
<p>Additional conditions can be handled using <code>elif</code> and <code>else</code>:</p>
</div>


<div class="in">
<pre>temp = data[&#39;temperature&#39;][0]
if temp &lt; 87:
    print &#34;The temperature is &lt; 87&#34;
elif temp &gt; 87:
    print &#34;The temperature is &gt; 87&#34;
else:
    print &#34; The temperature is equal to 87&#34;</pre>
</div>

<div class="out">
<pre> The temperature is equal to 87
</pre>
</div>

### Challenge


<div>
<p>Import the data from <code>mosquito_data_A2.csv</code>, determine the mean temperate, and loop over the temperature values. For each value print out whether it is greater than the mean, less than the mean, or equal to the mean.</p>
</div>

## Plotting


<div>
<p>The mathematician Richard Hamming once said, &quot;The purpose of computing is insight, not numbers,&quot; and the best way to develop insight is often to visualize data. The main plotting library in Python is <code>matplotlib</code>. To get started, let's tell the IPython Notebook that we want our plots displayed inline, rather than in a separate viewing window:</p>
</div>


<div class="in">
<pre>%matplotlib inline</pre>
</div>


<div>
<p>The <code>%</code> at the start of the line signals that this is a command for the notebook, rather than a statement in Python. Next, we will import the <code>pyplot</code> module from <code>matplotlib</code>, but since <code>pyplot</code> is a fairly long name to type repeatedly let's give it an alias.</p>
</div>


<div class="in">
<pre>from matplotlib import pyplot as plt</pre>
</div>


<div>
<p>This import statement shows two new things. First, we can import part of a library by using the <code>from library import submodule</code> syntax. Second, we can use a different name to refer to the imported library by using <code>as newname</code>.</p>
<p>Now, let's make a simple plot showing how the number of mosquitos varies over time. We'll use the site you've been doing exercises with since it has a longer time-series.</p>
</div>


<div class="in">
<pre>data = pandas.read_csv(&#39;mosquito_data_A2.csv&#39;)
plt.plot(data[&#39;year&#39;], data[&#39;mosquitos&#39;])</pre>
</div>

<div class="out">
<pre>[&lt;matplotlib.lines.Line2D at 0x4a88590&gt;]
<img src="../../intermediate/python/01-intro-python_files/intermediate/python/01-intro-python_66_1.png">
</pre>
</div>


<div>
<p>More complicated plots can be created by adding a little additional information. Let's say we want to look at how the different weather variables vary over time.</p>
</div>


<div class="in">
<pre>plt.figure(figsize=(10.0, 3.0))

plt.subplot(1, 2, 1)
plt.plot(data[&#39;year&#39;], data[&#39;temperature&#39;], &#39;ro-&#39;)
plt.xlabel(&#39;Year&#39;)
plt.ylabel(&#39;Temperature&#39;)

plt.subplot(1, 2, 2)
plt.plot(data[&#39;year&#39;], data[&#39;rainfall&#39;], &#39;bs-&#39;)
plt.xlabel(&#39;Year&#39;)
plt.ylabel(&#39;Rain Fall&#39;)
plt.show()</pre>
</div>

<div class="out">
<pre>
<img src="../../intermediate/python/01-intro-python_files/intermediate/python/01-intro-python_68_0.png">
</pre>
</div>

### Challenge


<div>
<p>Using the data in <code>mosquito_data_A2.csv</code> plot the relationship between the number of mosquitos and temperature and the number of mosquitos and rainfall.</p>
</div>

### Key Points


<div>
<ul>
<li>Import a library into a program using <code>import libraryname</code>.</li>
<li>Use the <code>pandas</code> library to work with data tables in Python.</li>
<li>Use <code>variable = value</code> to assign a value to a variable.</li>
<li>Use <code>print something</code> to display the value of <code>something</code>.</li>
<li>Use <code>dataframe['columnname']</code> to select a column of data.</li>
<li>Use <code>dataframe[start_row:stop_row]</code> to select rows from a data frame.</li>
<li>Indices start at 0, not 1.</li>
<li>Use <code>dataframe.mean()</code>, <code>dataframe.max()</code>, and <code>dataframe.min()</code> to calculate simple statistics.</li>
<li>Use <code>for x in list:</code> to loop over values</li>
<li>Use <code>if condition:</code> to make conditional decisions</li>
<li>Use the <code>pyplot</code> library from <code>matplotlib</code> for creating simple visualizations.</li>
</ul>
</div>

## Next steps


<div>
<p>With the requisite Python background out of the way, now we're ready to dig in to analyzing our data, and along the way learn how to write better code, more efficiently, that is more likely to be correct.</p>
</div>
