---
layout: lesson
root: ../..
---

# Modularization and Documentation


<div>
<p>Now that we've covered some of the basic syntax and libraries in Python we can start to tackle our data analysis problem. We are interested in understanding the relationship between the weather and the number of mosquitos so that we can plan mosquito control measures. Since we want to apply these mosquito control measures at a number of different sites we need to understand how the relationship varies across sites. Remember that we have a series of CSV files with each file containing the data for a single location.</p>
</div>

## Objectives


<div>
<ul>
<li>Write code for people, not computers</li>
<li>Break a program into chunks</li>
<li>Write and use functions in Python</li>
<li>Write useful documentation</li>
</ul>
</div>

## Starting small


<div>
<p>When approaching computational tasks like this one it is typically best to start small, check each piece of code as you go, and make incremental changes. This helps avoid marathon debugging sessions because it's much easier to debug one small piece of the code at a time than to write 100 lines of code and then try to figure out all of the different bugs in it.</p>
<p>Let's start by reading in the data from a single file and conducting a simple regression analysis on it. In fact, I would actually start by just importing the data and making sure that everything is coming in OK.</p>
</div>


<div class="in">
<pre>import pandas as pd

d = pd.read_csv(&#39;A2_mosquito_data.csv&#39;)
print d</pre>
</div>

<div class="out">
<pre>    year  temperature  rainfall  mosquitos
0   1960           82       200        180
1   1961           70       227        194
2   1962           89       231        207
3   1963           74       114        121
4   1964           78       147        140
5   1965           85       151        148
6   1966           86       172        162
7   1967           75       106        112
8   1968           70       276        230
9   1969           86       165        162
10  1970           83       222        198
11  1971           78       297        247
12  1972           87       288        248
13  1973           76       286        239
14  1974           86       231        202
15  1975           90       284        243
16  1976           76       190        175
17  1977           87       257        225
18  1978           88       128        133
19  1979           87       218        199
20  1980           81       206        184
21  1981           74       175        160
22  1982           85       202        187
23  1983           71       130        126
24  1984           80       225        200
25  1985           72       196        173
26  1986           76       261        222
27  1987           85       111        121
28  1988           83       247        210
29  1989           86       137        142
30  1990           82       159        152
31  1991           77       172        160
32  1992           74       280        231
33  1993           70       291        238
34  1994           77       126        125
35  1995           89       191        178
36  1996           83       298        248
37  1997           80       282        237
38  1998           86       219        195
39  1999           72       143        134
40  2000           79       262        221
41  2001           85       189        175
42  2002           86       205        186
43  2003           72       195        173
44  2004           78       148        146
45  2005           71       262        219
46  2006           88       255        226
47  2007           79       262        221
48  2008           73       198        176
49  2009           86       215        187
50  2010           87       127        129
</pre>
</div>


<div>
<p>The import seems to be working properly, so that's good news, but does anyone have anyone see anything about the code that they don't like?</p>
<p>That's right. The variable name I've chosen for the data doesn't really communicate any information to anyone about what it's holding, which means that when I come back to my code next month to change something I'm going to have a more difficult time understanding what the code is actually doing. This brings us to one of our first major lessons for the morning, which is that in order to understand what our code is doing so that we can quickly make changes in the future, we need to <em>write code for people, not computers</em>, and an important first step is to <em>use meaningful varible names</em>.</p>
</div>


<div class="in">
<pre>import pandas as pd

data = pd.read_csv(&#39;A2_mosquito_data.csv&#39;)
print data.head()</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
0  1960           82       200        180
1  1961           70       227        194
2  1962           89       231        207
3  1963           74       114        121
4  1964           78       147        140
</pre>
</div>


<div>
<p>The <code>.head()</code> method lets us just look at the first few rows of the data. A method is a function attached to an object that operates on that object. So in this case we can think of it as being equivalent to <code>head(data)</code>.</p>
<p>Everything looks good, but either global warming has gotten <em>really</em> out of control or the temperatures are in degrees Fahrenheit. Let's convert them to Celcius before we get started.</p>
<p>We don't need to reimport the data in our new cell because all of the executed cells in IPython Notebook share the same workspace. However, it's worth noting that if we close the notebook and then open it again it is necessary to rerun all of the individual blocks of code that a code block relies on before continuing. To rerun all of the cells in a notebook you can select <code>Cell -&gt; Run All</code> from the menu.</p>
</div>


<div class="in">
<pre>data[&#39;temperature&#39;] = (data[&#39;temperature&#39;] - 32) * 5 / 9.0
print data.head()</pre>
</div>

<div class="out">
<pre>   year  temperature  rainfall  mosquitos
0  1960    27.777778       200        180
1  1961    21.111111       227        194
2  1962    31.666667       231        207
3  1963    23.333333       114        121
4  1964    25.555556       147        140
</pre>
</div>


<div>
<p>That's better. Now let's go ahead and conduct a regression on the data. We'll use the <code>statsmodels</code> library to conduct the regression.</p>
</div>


<div class="in">
<pre>import statsmodels.api as sm

regr_results = sm.OLS.from_formula(&#39;mosquitos ~ temperature + rainfall&#39;, data).fit()
regr_results.summary()</pre>
</div>

<div class="out">
<pre>&lt;class &#39;statsmodels.iolib.summary.Summary&#39;&gt;
&#34;&#34;&#34;
                            OLS Regression Results                            
==============================================================================
Dep. Variable:              mosquitos   R-squared:                       0.997
Model:                            OLS   Adj. R-squared:                  0.997
Method:                 Least Squares   F-statistic:                     7889.
Date:                Sat, 18 Jan 2014   Prob (F-statistic):           3.68e-61
Time:                        13:10:18   Log-Likelihood:                -111.54
No. Observations:                  51   AIC:                             229.1
Df Residuals:                      48   BIC:                             234.9
Df Model:                           2                                         
===============================================================================
                  coef    std err          t      P&gt;|t|      [95.0% Conf. Int.]
-------------------------------------------------------------------------------
Intercept      17.5457      2.767      6.341      0.000        11.983    23.109
temperature     0.8719      0.092      9.457      0.000         0.687     1.057
rainfall        0.6967      0.006    125.385      0.000         0.686     0.708
==============================================================================
Omnibus:                        1.651   Durbin-Watson:                   1.872
Prob(Omnibus):                  0.438   Jarque-Bera (JB):                0.906
Skew:                          -0.278   Prob(JB):                        0.636
Kurtosis:                       3.343   Cond. No.                     1.92e+03
==============================================================================

Warnings:
[1] The condition number is large, 1.92e+03. This might indicate that there are
strong multicollinearity or other numerical problems.
&#34;&#34;&#34;</pre>
</div>


<div>
<p>As you can see <code>statsmodels</code> lets us use the names of the columns in our dataframe to clearly specify the form of the statistical model we want to fit. This also makes the code more readable since the model we are fitting is written in a nice, human readable, manner. The <code>summary</code> method gives us a visual representation of the results. This summary is nice to look at, but it isn't really useful for doing more computation, so we can look up particular values related to the regression using the <code>regr_results</code> attributes. These are variables that are attached to <code>regr_results</code>.</p>
</div>


<div class="in">
<pre>print regr_results.params
print regr_results.rsquared</pre>
</div>

<div class="out">
<pre>Intercept      17.545739
temperature     0.871943
rainfall        0.696717
dtype: float64
0.996966873691
</pre>
</div>


<div>
<p>If we want to hold onto these values for later we can assign them to variables:</p>
</div>


<div class="in">
<pre>parameters = regr_results.params
rsquared = regr_results.rsquared</pre>
</div>


<div>
<p>And then we can plot the observed data against the values predicted by our regression to visualize the results. First, remember to tell the notebook that we want our plots to appear in the notebook itself.</p>
</div>


<div class="in">
<pre>%matplotlib inline</pre>
</div>


<div class="in">
<pre>import matplotlib.pyplot as plt

predicted = parameters[0] + parameters[1] * data[&#39;temperature&#39;] + parameters[2] * data[&#39;rainfall&#39;]
plt.plot(predicted, data[&#39;mosquitos&#39;], &#39;ro&#39;)
min_mosquitos, max_mosquitos = min(data[&#39;mosquitos&#39;]), max(data[&#39;mosquitos&#39;])
plt.plot([min_mosquitos, max_mosquitos], [min_mosquitos, max_mosquitos], &#39;k-&#39;)</pre>
</div>

<div class="out">
<pre>[&lt;matplotlib.lines.Line2D at 0x56eb950&gt;]
<img src="../../intermediate/python/02-modularization-documentation_files/intermediate/python/02-modularization-documentation_19_1.png">
</pre>
</div>


<div>
<p>OK, great. So putting this all together we now have a piece of code that imports the modules we need, loads the data into memory, fits a regression to the data, and stores the parameters and fit of data.</p>
</div>


<div class="in">
<pre>import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

data = pd.read_csv(&#39;A2_mosquito_data.csv&#39;)
data[&#39;temperature&#39;] = (data[&#39;temperature&#39;] - 32) * 5 / 9.0
regr_results = sm.OLS.from_formula(&#39;mosquitos ~ temperature + rainfall&#39;, data).fit()
parameters = regr_results.params
rsquared = regr_results.rsquared
predicted = parameters[0] + parameters[1] * data[&#39;temperature&#39;] + parameters[2] * data[&#39;rainfall&#39;]
plt.plot(predicted, data[&#39;mosquitos&#39;], &#39;ro&#39;)
min_mosquitos, max_mosquitos = min(data[&#39;mosquitos&#39;]), max(data[&#39;mosquitos&#39;])
plt.plot([min_mosquitos, max_mosquitos], [min_mosquitos, max_mosquitos], &#39;k-&#39;)
print parameters
print &#34;R^2 = &#34;, rsquared</pre>
</div>

<div class="out">
<pre>Intercept      17.545739
temperature     0.871943
rainfall        0.696717
dtype: float64
R^2 =  0.996966873691

<img src="../../intermediate/python/02-modularization-documentation_files/intermediate/python/02-modularization-documentation_21_1.png">
</pre>
</div>

## Functions


<div>
<p>The next thing we need to do is loop over all of the possible data files, but in order to do that we're going to need to grow our code some more. Since our brain can only easily hold 5-7 pieces of information at once, and our code already has more than that many pieces, we need to start breaking our code into manageable sized chunks. This will let us read and understand the code more easily and make it easier to reuse pieces of our code. We'll do this using functions.</p>
<p>Functions in Python take the general form</p>
<pre class="sourceCode python"><code class="sourceCode python"><span class="kw">def</span> function_name(inputs):
    do stuff
    <span class="kw">return</span> output</code></pre>
<p>So, if we want to write a function that returns the value of a number squared we could use:</p>
</div>


<div class="in">
<pre>def square(x):
    x_squared = x ** 2
    return x_squared

print &#34;Four squared is&#34;, square(4)
print &#34;Five squared is&#34;, square(5)</pre>
</div>

<div class="out">
<pre>Four squared is 16
Five squared is 25
</pre>
</div>


<div>
<p>We can also just return the desired value directly.</p>
</div>


<div class="in">
<pre>def square(x):
    return x ** 2

print square(3)</pre>
</div>

<div class="out">
<pre>9
</pre>
</div>


<div>
<p>And remember, if we want to use the result of the function later we need to store it somewhere.</p>
</div>


<div class="in">
<pre>two_squared = square(2)
print two_squared</pre>
</div>

<div class="out">
<pre>4
</pre>
</div>

### Challenges


<div>
<p>1. Write a function that converts temperature from Fahrenheit to Celcius and use it to replace</p>
<pre class="sourceCode python"><code class="sourceCode python">data[<span class="st">&#39;temperature&#39;</span>] = (data[<span class="st">&#39;temperature&#39;</span>] - <span class="dv">32</span>) * <span class="dv">5</span> / <span class="fl">9.0</span></code></pre>
<p>in our program.</p>
<p>2. Write a function called <code>analyze()</code> that takes <code>data</code> as an input, performs the regression, makes the observed-predicted plot, and returns <code>parameters</code>.</p>
<p><em>Walk through someone's result. When discussing talk about different names. E.g., fahr_to_celcius is better than temp_to_celcius since it is explicit both the input and the output. Talk about the fact that even though this doesn't save us any lines of code it's still easier to read.</em></p>
</div>

## The call stack


<div>
<p>Let's take a closer look at what happens when we call fahr_to_celsius(32.0). To make things clearer, we'll start by putting the initial value 32.0 in a variable and store the final result in one as well:</p>
</div>


<div class="in">
<pre>def fahr_to_celsius(tempF):
    tempC = (tempF - 32) * 5 / 9.0
    return tempC

original = 32.0
final = fahr_to_celsius(original)</pre>
</div>


<div>
<h4 id="call-stack-initial-state">Call Stack (Initial State)</h4>
<p>When the first three lines of this function are executed the function is created, but nothing happens. The function is like a recipe, it contains the information about how to do something, but it doesn't do so until you explicitly ask it to. We then create the variable <code>original</code> and assign the value 32.0 to it. The values <code>tempF</code> and <code>tempC</code> don't currently exist.</p>
<h4 id="call-stack-immediately-after-function-call">Call Stack Immediately After Function Call</h4>
<p>When we call <code>fahr_to_celsius</code>, Python creates another stack frame to hold fahr_to_celsius's variables. Upon creation this stack frame only includes the inputs being passed to the function, so in our case <code>tempF</code>. As the function is executed variables created by the function are stored in the functions stack frame, so <code>tempC</code> is created in the <code>fahr_to_celsius</code> stack frame.</p>
<h4 id="call-stack-at-end-of-function-call">Call Stack At End Of Function Call</h4>
<p>When the call to <code>fahr_to_celsius</code> returns a value, Python throws away <code>fahr_to_celsius</code>'s stack frame, including all of the variables it contains, and creates a new variable in the original stack frame to hold the temperature in Celsius.</p>
<h4 id="call-stack-after-end">Call Stack After End</h4>
<p>This final stack frame is always there; it holds the variables we defined outside the functions in our code. What it doesn't hold is the variables that were in the other stack frames. If we try to get the value of <code>tempF</code> or <code>tempC</code> after our functions have finished running, Python tells us that there's no such thing:</p>
</div>


<div class="in">
<pre>print tempC</pre>
</div>

<div class="out">
<pre>---------------------------------------------------------------------------
NameError                                 Traceback (most recent call last)
&lt;ipython-input-14-3054d7679e45&gt; in &lt;module&gt;()
----&gt; 1 print tempC

NameError: name &#39;tempC&#39; is not defined</pre>
</div>


<div>
<p>The reason for this is encapsulation, and it's one of the key to writing correct, comprehensible programs. A function's job is to turn several operations into one so that we can think about a single function call instead of a dozen or a hundred statements each time we want to do something. That only works if functions don't interfere with each other by potentially changing the same variables; if they do, we have to pay attention to the details once again, which quickly overloads our short-term memory.</p>
</div>

## Testing Functions


<div class="">
<p>Once we start putting things into functions so that we can re-use them, we need to start testing that those functions are working correctly. The most basic thing we can do is some informal testing to make sure the function is doing what it is supposed to do. To see how to do this, let's write a function to center the values in a dataset prior to conducting statistical analysis. Centering means setting the mean of each variable to be the same value, typically zero.</p>
</div>


<div class="in">
<pre>def center(data):
    return data - data.mean()</pre>
</div>


<div class="">
<p>We could test this on our actual data, but since we don't know what the values ought to be, it will be hard to tell if the result was correct. Instead, let's create a made up data frame where we know what the result should look like.</p>
</div>


<div class="in">
<pre>import pandas as pd

test_data = pd.DataFrame([[1, 1], [1, 2]])
print test_data</pre>
</div>

<div class="out">
<pre>   0  1
0  1  1
1  1  2
</pre>
</div>


<div>
<p>Now that we've made some test data we need to figure out what we think the result should be and we need to do this <em>before</em> we run the test. This is important because we are biased to believe that any result we get back is correct, and we want to avoid that bias. This also helps make sure that we are confident in what we want the code to do. So, what should the result of running <code>center(data)</code> be?</p>
<p>OK, let's go ahead and run the function.</p>
</div>


<div class="in">
<pre>print center(test_data)</pre>
</div>

<div class="out">
<pre>   0    1
0  0 -0.5
1  0  0.5
</pre>
</div>


<div class="">
<p>That looks right, so let's try <code>center</code> on our real data:</p>
</div>


<div class="in">
<pre>data = pd.read_csv(&#39;A2_mosquito_data.csv&#39;)
print center(data)</pre>
</div>

<div class="out">
<pre>    year  temperature    rainfall  mosquitos
0    -25     1.607843   -7.039216  -5.235294
1    -24   -10.392157   19.960784   8.764706
2    -23     8.607843   23.960784  21.764706
3    -22    -6.392157  -93.039216 -64.235294
4    -21    -2.392157  -60.039216 -45.235294
5    -20     4.607843  -56.039216 -37.235294
6    -19     5.607843  -35.039216 -23.235294
7    -18    -5.392157 -101.039216 -73.235294
8    -17   -10.392157   68.960784  44.764706
9    -16     5.607843  -42.039216 -23.235294
10   -15     2.607843   14.960784  12.764706
11   -14    -2.392157   89.960784  61.764706
12   -13     6.607843   80.960784  62.764706
13   -12    -4.392157   78.960784  53.764706
14   -11     5.607843   23.960784  16.764706
15   -10     9.607843   76.960784  57.764706
16    -9    -4.392157  -17.039216 -10.235294
17    -8     6.607843   49.960784  39.764706
18    -7     7.607843  -79.039216 -52.235294
19    -6     6.607843   10.960784  13.764706
20    -5     0.607843   -1.039216  -1.235294
21    -4    -6.392157  -32.039216 -25.235294
22    -3     4.607843   -5.039216   1.764706
23    -2    -9.392157  -77.039216 -59.235294
24    -1    -0.392157   17.960784  14.764706
25     0    -8.392157  -11.039216 -12.235294
26     1    -4.392157   53.960784  36.764706
27     2     4.607843  -96.039216 -64.235294
28     3     2.607843   39.960784  24.764706
29     4     5.607843  -70.039216 -43.235294
30     5     1.607843  -48.039216 -33.235294
31     6    -3.392157  -35.039216 -25.235294
32     7    -6.392157   72.960784  45.764706
33     8   -10.392157   83.960784  52.764706
34     9    -3.392157  -81.039216 -60.235294
35    10     8.607843  -16.039216  -7.235294
36    11     2.607843   90.960784  62.764706
37    12    -0.392157   74.960784  51.764706
38    13     5.607843   11.960784   9.764706
39    14    -8.392157  -64.039216 -51.235294
40    15    -1.392157   54.960784  35.764706
41    16     4.607843  -18.039216 -10.235294
42    17     5.607843   -2.039216   0.764706
43    18    -8.392157  -12.039216 -12.235294
44    19    -2.392157  -59.039216 -39.235294
45    20    -9.392157   54.960784  33.764706
46    21     7.607843   47.960784  40.764706
47    22    -1.392157   54.960784  35.764706
48    23    -7.392157   -9.039216  -9.235294
49    24     5.607843    7.960784   1.764706
50    25     6.607843  -80.039216 -56.235294
</pre>
</div>


<div class="">
<p>It's hard to tell from the default output whether the result is correct, but there are a few simple tests that will reassure us:</p>
</div>


<div class="in">
<pre>print &#39;original mean:&#39;
print data.mean()
centered = center(data)
print
print &#39;mean of centered data:&#39;
centered.mean()</pre>
</div>

<div class="out">
<pre>original mean:
year           1985.000000
temperature      80.392157
rainfall        207.039216
mosquitos       185.235294
dtype: float64

mean of centered data:
year           0.000000e+00
temperature    1.393221e-15
rainfall       6.687461e-15
mosquitos     -1.337492e-14
dtype: float64</pre>
</div>


<div class="">
<p>The mean of the centered data is very close to zero; it's not quite zero because of floating point precision issues. We can even go further and check that the standard deviation hasn't changed (which it shouldn't if we've just centered the data):</p>
</div>


<div class="in">
<pre>print &#39;std dev before and after:&#39;
print data.std()
print
print centered.std()</pre>
</div>

<div class="out">
<pre>std dev before and after:
year           14.866069
temperature     6.135400
rainfall       56.560396
mosquitos      39.531551
dtype: float64

year           14.866069
temperature     6.135400
rainfall       56.560396
mosquitos      39.531551
dtype: float64
</pre>
</div>


<div class="">
<p>The standard deviations look the same. It's still possible that our function is wrong, but it seems unlikely enough that we we're probably in good shape for now.</p>
</div>

## Documentation


<div>
<p>OK, the <code>center</code> function seems to be working fine. Does anyone else see anything that's missing before we move on?</p>
<p>Yes, we should write some <a href="../../gloss.html#documentation">documentation</a> to remind ourselves later what it's for and how to use it. This function may be fairly straightforward, but in most cases it won't be so easy to remember exactly what a function is doing in a few months. Just imagine looking at our <code>analyze</code> function a few months in the future and trying to remember exactly what it was doing just based on the code.</p>
<p>The usual way to put documentation in code is to add <a href="../../gloss.html#comment">comments</a> like this:</p>
</div>


<div class="in">
<pre># center(data): return a new DataFrame containing the original data centered around zero.
def center(data, desired):
    return data - data.mean()</pre>
</div>


<div class="">
<p>There's a better way to do this in Python. If the first thing in a function is a string that isn't assigned to a variable, that string is attached to the function as its documentation:</p>
</div>


<div class="in">
<pre>def center(data, desired):
    &#34;&#34;&#34;Return a new DataFrame containing the original data centered around zero.&#34;&#34;&#34;
    return data - data.mean()</pre>
</div>


<div class="">
<p>This is better because we can now ask Python's built-in help system to show us the documentation for the function.</p>
</div>


<div class="in">
<pre>help(center)</pre>
</div>

<div class="out">
<pre>Help on function center in module __main__:

center(data, desired)
    Return a new DataFrame containing the original data centered around zero.

</pre>
</div>


<div class="">
<p>A string like this is called a <a href="../../gloss.html#docstring">docstring</a> and there are also automatic documentation generators that use these docstrings to produce documentation for users. We use triple quotes because it allows us to include multiple lines of text and because it is considered good Python style.</p>
</div>


<div class="in">
<pre>def center(data):
    &#34;&#34;&#34;Return a new array containing the original data centered on zero
    
    Example:
    &gt;&gt;&gt; import pandas
    &gt;&gt;&gt; data = pandas.DataFrame([[0, 1], [0, 2])
    &gt;&gt;&gt; center(data)
       0    1
    0  0 -0.5
    1  0  0.5

     
    &#34;&#34;&#34;
    return data - data.mean()

help(center)</pre>
</div>

<div class="out">
<pre>Help on function center in module __main__:

center(data)
    Return a new array containing the original data centered on zero
    
    Example:
    &gt;&gt;&gt; import pandas
    &gt;&gt;&gt; data = pandas.DataFrame([[0, 1], [0, 2])
    &gt;&gt;&gt; center(data)
       0    1
    0  0 -0.5
    1  0  0.5

</pre>
</div>

### Challenge


<div>
<ol style="list-style-type: decimal">
<li>Test your temperature conversion function to make sure it's working (think about some temperatures that you easily know the conversion for).</li>
<li>Add documentation to both the temperature conversation function and the analysis function.</li>
</ol>
</div>

## Looping over files


<div>
<p>So now our code looks something like this:</p>
</div>


<div class="in">
<pre>import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

def fahr_to_celsius(tempF):
    &#34;&#34;&#34;Convert fahrenheit to celsius&#34;&#34;&#34;
    tempC = (tempF - 32) * 5 / 9.0
    return tempC

def analyze(data):
    &#34;&#34;&#34;Perform regression analysis on mosquito data
    
    Takes a dataframe as input that includes columns named &#39;temperature&#39;,
    &#39;rainfall&#39;, and &#39;mosquitos&#39;.
    
    Performs a multiple regression to predict the number of mosquitos.
    Creates an observed-predicted plot of the result and
    returns the parameters of the regression.
    
    &#34;&#34;&#34;
    regr_results = sm.OLS.from_formula(&#39;mosquitos ~ temperature + rainfall&#39;, data).fit()
    parameters = regr_results.params
    predicted = parameters[0] + parameters[1] * data[&#39;temperature&#39;] + parameters[2] * data[&#39;rainfall&#39;]
    plt.figure()
    plt.plot(predicted, data[&#39;mosquitos&#39;], &#39;ro&#39;)
    min_mosquitos, max_mosquitos = min(data[&#39;mosquitos&#39;]), max(data[&#39;mosquitos&#39;])
    plt.plot([min_mosquitos, max_mosquitos], [min_mosquitos, max_mosquitos], &#39;k-&#39;)
    return parameters

data = pd.read_csv(&#39;A2_mosquito_data.csv&#39;)
data[&#39;temperature&#39;] = fahr_to_celsius(data[&#39;temperature&#39;])
regr_results = analyze(data)
print parameters</pre>
</div>

<div class="out">
<pre>Intercept      17.545739
temperature     0.871943
rainfall        0.696717
dtype: float64

<img src="../../intermediate/python/02-modularization-documentation_files/intermediate/python/02-modularization-documentation_64_1.png">
</pre>
</div>


<div>
<p>Now we want to loop over all of the possible data files, and to do that we need to know their names. If we only had a dozen files we could write them all down, but if we have hundreds of files or the filenames change then that won't really work. Fortunately Python has a built in library to help us find the files we want to work with called <code>glob</code>.</p>
</div>


<div class="in">
<pre>import glob

filenames = glob.glob(&#39;*.csv&#39;)
print filenames</pre>
</div>

<div class="out">
<pre>[&#39;A1_mosquito_data.csv&#39;, &#39;B2_mosquito_data.csv&#39;, &#39;A3_mosquito_data.csv&#39;, &#39;B1_mosquito_data.csv&#39;, &#39;A2_mosquito_data.csv&#39;]
</pre>
</div>


<div>
<p>The object returned by <code>glob</code> is a list of strings. A list is a Python data type that holds a group of potentially heterogenous values. That means it can hold pretty much anything, including functions.</p>
</div>


<div class="in">
<pre>mylist = [1, &#39;a&#39;, center]
print mylist</pre>
</div>

<div class="out">
<pre>[1, &#39;a&#39;, &lt;function center at 0x5c63b18&gt;]
</pre>
</div>


<div>
<p>In this case all of the values are strings that contain the names of all of the files that match the expression given to <code>glob</code>, so in this case all of the files with the <code>.csv</code> extension.</p>
<p>Let's restrict the filenames a little more finely, so that we don't accidentally get any data we don't want, and print out the filenames one at a time.</p>
</div>


<div class="in">
<pre>filenames =glob.glob(&#39;*data.csv&#39;)
for filename in filenames:
    print filename</pre>
</div>

<div class="out">
<pre>A1_mosquito_data.csv
B2_mosquito_data.csv
A3_mosquito_data.csv
B1_mosquito_data.csv
A2_mosquito_data.csv
</pre>
</div>

### Challenge


<div>
<p>Modify your code to loop over all of the files in your directory, making an observed-predicted plot for each file and printing the parameters.</p>
</div>
