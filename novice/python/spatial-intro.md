---
layout: lesson
root: ../..
---

## Visualizing Spatial Data


<div>
<blockquote>
<h2>A note to students and instructors</h2>
<p>This lesson requires the <a href="http://matplotlib.org/basemap">Basemap</a> toolkit for Matplotlib. This library is not distributed with Matplotlib; if you are using Continuum's Anaconda distribution, you can obtain it by running:</p>
<pre><code>conda install basemap</code></pre>
<p>If you are using Enthought Canopy and have the full version or an academic license, Basemap should already be installed on your system. Otherwise, you will need to follow the <a href="http://matplotlib.org/basemap/users/installing.html">installation instructions</a> in the Basemap documentation. Using one of the two scientific distributions is preferred in most instances.</p>
</blockquote>
</div>

## Visualizing spatial data


<div>
<p>We are examining some simple spatial coordinate data, specifically the location of all of the previous Software Carpentry bootcamps. The data set is stored in <a href="../../gloss.html#csv">comma-separated values</a> (CSV) format. After the header line marked with a <code>#</code>, each row contains the latitude and longitude for each bootcamp, separated by a comma:</p>
<pre><code># Latitude, Longitude
43.661476,-79.395189
39.332604,-76.623190
45.703255, 13.718013
43.661476,-79.395189
39.166381,-86.526621
...</code></pre>
<p>We want to:</p>
<ul>
<li>load the data into our analysis environment;</li>
<li>inspect the data; and</li>
<li>visualize it.</li>
</ul>
<p>To do this, we'll delve into working with Python.</p>
</div>


<div>
<h2 id="objectives">Objectives</h2>
<ul>
<li>Explain what a library is, and what libraries are used for.</li>
<li>Load a Python library and use the things it contains.</li>
<li>Read tabular data from a file into a program.</li>
<li>Display simple visualizations of the data</li>
</ul>
</div>

### Loading the Data


<div>
<p>In order to work with the coordinates stored in the file, we need to <a href="../../gloss.html#import">import</a> a library called NumPy that is designed to handle arrays of data.</p>
</div>


<div class="in">
<pre>import numpy as np</pre>
</div>


<div>
<p>It's very common to create an <a href="../../gloss.html#alias-library">alias</a> for a library when importing it in order to reduce the amount of typing we have to do. We can now refer to this library in the code as <code>np</code> instead of typing out <code>numpy</code> each time we want to use it.</p>
<p>We can now ask NumPy to read our data file:</p>
</div>


<div class="in">
<pre>lat, lon = np.loadtxt(&#39;swc_bc_coords.csv&#39;, delimiter=&#39;,&#39;, unpack=True)</pre>
</div>


<div>
<p>The expression <code>np.loadtxt(...)</code> means, &quot;Run the function <code>loadtxt</code> that belongs to the <code>numpy</code> library.&quot; This <a href="../../gloss.html#dotted-notation">dotted notation</a> is used everywhere in Python to refer to the parts of larger things.</p>
<p><code>np.loadtxt</code> has three <a href="../../gloss.html#parameter">parameters</a>:</p>
<ol style="list-style-type: decimal">
<li>the name of the file we want to read,</li>
<li>the <a href="../../gloss.html#delimiter">delimiter</a> that separates values on a line, and</li>
<li>a <a href="../../gloss.html#flag">flag</a> called <code>unpack</code>.</li>
</ol>
<p>The first two parameters both need to be <a href="../../gloss.html#string">character strings</a> (or strings for short), so we put them in quotes. Setting <code>unpack</code> to <code>True</code> tells <code>np.loadtxt</code> to take the first and second column of data and give them back to us separately so that we can <a href="../../gloss.html#assignment">assign</a> them to the two <a href="../../gloss.html#variable">variables</a> <code>lat</code> and <code>lon</code>.</p>
<blockquote>
<h4>Trying to be Helpful</h4>
<p><code>np.loadtxt</code> automatically skips the line with the header information (the one starting with '#'), since it recognizes that this line is a <a href="../../gloss.html#comment">comment</a> and does not contain numerical data.</p>
</blockquote>
<p>When we are finished typing and press Shift+Enter, the notebook runs our command. <code>lat</code> and <code>lon</code> now contain our data, which we can inspect by printing either of the variables:</p>
</div>


<div class="in">
<pre>print lat</pre>
</div>

<div class="out">
<pre>[ 43.661476  39.332604  45.703255  43.661476  39.166381  36.802151
  37.808381  41.790113  41.744949  51.559882  42.727288  54.980095
  53.523454  49.261715  39.32758   48.831673  42.359133  43.47013
  44.632261  43.783551  53.948193  59.939959  40.808078  40.428267
  37.875928  49.261715  37.8695    54.980095  34.141411  38.831513
  51.757137  43.261328  38.648056  32.89533   34.227425  21.300662
  55.945328  30.283599  49.261715  41.790113  45.417417  43.469128
  49.261715  48.264934  43.647118  48.53698   40.808078  37.228384
  49.261715 -33.773636 -37.825328  47.655965  37.875928  38.031441
  33.900058  41.744949  22.3101    32.236358  51.524789 -33.929492
  53.467102  37.8695    53.478349  48.82629   39.291389  43.07718   52.33399
  54.32707   39.07141   37.42949   37.875928  43.64712   51.759865
  38.54926   36.00803   50.060833  36.00283   40.03131   42.388889
  53.52345   50.937716  42.35076   41.789722  49.276765  32.887151
  41.790113  42.3625    30.283599 -43.523333  35.20859   59.939959
  30.538978  39.166381  51.377743  37.228384  41.7408    41.70522   47.655
  40.443322  44.968657  38.958455  32.30192   43.07718   41.66293
  51.457971  43.468889  42.724085 -34.919159  49.261111 -37.9083    34.052778
  41.526667]
</pre>
</div>


<div>
<h2 id="visualizing-the-data">Visualizing the Data</h2>
</div>


<div>
<p>A lot of tools for working with data are built into NumPy's arrays. We'll explore some of them in later lessons, but for now let's just make a simple plot of the data using another library called <code>matplotlib</code>. First, let's tell the IPython Notebook that we want our plots displayed inline, rather than in a separate viewing window:</p>
</div>


<div class="in">
<pre>%matplotlib inline</pre>
</div>


<div>
<p>The <code>%</code> at the start of the line signals that this is a command for the notebook, rather than a statement in Python. Next, we will import the <code>pyplot</code> module from <code>matplotlib</code> (again using an alias to cut down on typing) and use one of the commands it defines to make plot a point for each latitude, longitude pair of data.</p>
</div>


<div class="in">
<pre>from matplotlib import pyplot as plt
plt.plot(lon, lat, &#39;o&#39;)</pre>
</div>

<div class="out">
<pre>[&lt;matplotlib.lines.Line2D at 0x10622ce50&gt;]</pre>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li>Plot the dots with a different color according to the continent they would be on.</li>
</ol>
</div>

## Mapping the Data


<div>
<p>While matplotlib provides a simple facility for visualizing numerical data in a variety of ways, we will use a supplementary toolkit called <code>Basemap</code> that enhances matplotlib to specifically deal with spatial data. We need to import this library and can do so using:</p>
</div>


<div class="in">
<pre>from mpl_toolkits.basemap import Basemap</pre>
</div>


<div>
<p>Now let's create a Basemap object that will allow us to project the coordinates onto map. For this example we will use a <a href="http://en.wikipedia.org/wiki/Robinson_projection">Robinson Projection</a>:</p>
</div>


<div class="in">
<pre>map = Basemap(projection=&#39;robin&#39;, lat_0=0.0, lon_0=0.0)</pre>
</div>


<div>
<p>The <code>projection</code> parameter is self-explanatory; the parameters <code>lat_0</code> and <code>lon_0</code> define the center of the map.</p>
<p>Now that we have a map, we can use <code>pyplot</code> to create a figure to display things in, and then tell the map which features to display: coastlines, country borders, shaded continents, and so on.</p>
<p>FIXME: explain the <code>x, y = map(lon, lat)</code> call.</p>
</div>


<div class="in">
<pre>plt.figure(figsize=(12,12))
map.drawcoastlines()
map.drawcountries()
map.fillcontinents()
map.drawmeridians(np.arange(-180,180,20))
map.drawparallels(np.arange(-90,90,20))

x, y = map(lon, lat)
map.plot(x, y, &#39;o&#39;, markersize=4, color=&#39;red&#39;)</pre>
</div>

<div class="out">
<pre>[&lt;matplotlib.lines.Line2D at 0x108c88f10&gt;]</pre>
</div>


<div>
<p>The final line of the above code cell mimics matplotlib's built-in <code>plot</code> method to plot our projected coordinates onto the map.</p>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>Integrate the coloring scheme from Exercise 1 into the Basemap projection.</p></li>
<li><p>Try out a different projection that better shows the boot camp locations in North America. (You can find a list of projections in the <a href="http://matplotlib.org/basemap/users/mapsetup.html">Basemap documentation</a>.)</p></li>
</ol>
</div>
