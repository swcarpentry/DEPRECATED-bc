---
layout: lesson
root: ../..
---

## Calculating New Values


<div class="objectives">
<h4 id="objectives">Objectives</h4>
<ul>
<li>Write queries that calculate new values for each selected record.</li>
</ul>
</div>


<div>
<p>After carefully re-reading the expedition logs, we realize that the radiation measurements they report may need to be corrected upward by 5%. Rather than modifying the stored data, we can do this calculation on the fly as part of our query:</p>
</div>


<div class="in">
<pre>%load_ext sqlitemagic</pre>
</div>


<div class="in">
<pre>%%sqlite survey.db
select 1.05 * reading from Survey where quant=&#39;rad&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>10.311</td></tr>
<tr><td>8.19</td></tr>
<tr><td>8.8305</td></tr>
<tr><td>7.581</td></tr>
<tr><td>4.5675</td></tr>
<tr><td>2.2995</td></tr>
<tr><td>1.533</td></tr>
<tr><td>11.8125</td></tr>
</table></pre>
</div>


<div>
<p>When we run the query, the expression <code>1.05 * reading</code> is evaluated for each row. Expressions can use any of the fields, all of usual arithmetic operators, and a variety of common functions. (Exactly which ones depends on which database manager is being used.) For example, we can convert temperature readings from Fahrenheit to Celsius and round to two decimal places:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select taken, round(5*(reading-32)/9, 2) from Survey where quant=&#39;temp&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>734</td><td>-29.72</td></tr>
<tr><td>735</td><td>-32.22</td></tr>
<tr><td>751</td><td>-28.06</td></tr>
<tr><td>752</td><td>-26.67</td></tr>
</table></pre>
</div>


<div>
<p>We can also combine values from different fields, for example by using the string concatenation operator <code>||</code>:</p>
</div>


<div class="in">
<pre>%%sqlite survey.db
select personal || &#39; &#39; || family from Person;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>William Dyer</td></tr>
<tr><td>Frank Pabodie</td></tr>
<tr><td>Anderson Lake</td></tr>
<tr><td>Valentina Roerich</td></tr>
<tr><td>Frank Danforth</td></tr>
</table></pre>
</div>


<div>
<blockquote>
<p>It may seem strange to use <code>personal</code> and <code>family</code> as field names instead of <code>first</code> and <code>last</code>, but it's a necessary first step toward handling cultural differences. For example, consider the following rules:</p>
</blockquote>
<table>
  <tr> <th>
Full Name
</th> <th>
Alphabetized Under
</th> <th>
Reason
</th> </tr>
  <tr> <td>
Liu Xiaobo
</td> <td>
Liu
</td> <td>
Chinese family names come first
</td> </tr>
  <tr> <td> 
Leonardo da Vinci
</td> <td>
Leonardo
</td> <td>
&quot;da Vinci&quot; just means &quot;from Vinci&quot;
</td> </tr>
  <tr> <td> 
Catherine de Medici
</td> <td>
Medici
</td> <td>
family name
</td> </tr>
  <tr> <td> 
Jean de La Fontaine
</td> <td>
La Fontaine
</td> <td>
family name is &quot;La Fontaine&quot;
</td> </tr>
  <tr> <td> 
Juan Ponce de Leon
</td> <td>
Ponce de Leon
</td> <td>
full family name is &quot;Ponce de Leon&quot;
</td> </tr>
  <tr> <td> 
Gabriel Garcia Marquez
</td> <td>
Garcia Marquez
</td> <td>
double-barrelled Spanish surnames
</td> </tr>
  <tr> <td> 
Wernher von Braun
</td> <td>
von <em>or</em> Braun
</td> <td>
depending on whether he was in Germany or the US
</td> </tr>
  <tr> <td> 
Elizabeth Alexandra May Windsor
</td> <td>
Elizabeth
</td> <td>
monarchs alphabetize by the name under which they reigned
</td> </tr>
  <tr> <td> 
Thomas a Beckett
</td> <td>
Thomas
</td> <td>
and saints according to the names by which they were canonized
</td> </tr>
</table>

<blockquote>
<p>Clearly, even a two-part division into &quot;personal&quot; and &quot;family&quot; isn't enough...</p>
</blockquote>
</div>


<div>
<h4 id="challenges">Challenges</h4>
<ol style="list-style-type: decimal">
<li><p>After further reading, we realize that Valentina Roerich was reporting salinity as percentages. Write a query that returns all of her salinity measurements from the <code>Survey</code> table with the values divided by 100.</p></li>
<li><p>The <code>union</code> operator combines the results of two queries:</p></li>
</ol>
</div>


<div class="in">
<pre>%%sqlite survey.db
select * from Person where ident=&#39;dyer&#39; union select * from Person where ident=&#39;roe&#39;;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>dyer</td><td>William</td><td>Dyer</td></tr>
<tr><td>roe</td><td>Valentina</td><td>Roerich</td></tr>
</table></pre>
</div>


<div>
<p>Use <code>union</code> to create a consolidated list of salinity measurements in which Roerich's, and only Roerich's, have been corrected as described in the previous challenge. The output should be something like:</p>
<table>
  <tr> <td>
619
</td> <td>
0.13
</td> </tr>
  <tr> <td>
622
</td> <td>
0.09
</td> </tr>
  <tr> <td>
734
</td> <td>
0.05
</td> </tr>
  <tr> <td>
751
</td> <td>
0.1
</td> </tr>
  <tr> <td>
752
</td> <td>
0.09
</td> </tr>
  <tr> <td>
752
</td> <td>
0.416
</td> </tr>
  <tr> <td>
837
</td> <td>
0.21
</td> </tr>
  <tr> <td>
837
</td> <td>
0.225
</td> </tr>
</table>


</div>


<div>
<ol start="3" style="list-style-type: decimal">
<li>The site identifiers in the <code>Visited</code> table have two parts separated by a '-':</li>
</ol>
</div>


<div class="in">
<pre>%%sqlite survey.db
select distinct site from Visited;</pre>
</div>

<div class="out">
<pre><table>
<tr><td>DR-1</td></tr>
<tr><td>DR-3</td></tr>
<tr><td>MSK-4</td></tr>
</table></pre>
</div>


<div>
<p>Some major site identifiers are two letters long and some are three. The &quot;in string&quot; function <code>instr(X, Y)</code> returns the 1-based index of the first occurrence of string Y in string X, or 0 if Y does not exist in X. The substring function <code>substr(X, I)</code> returns the substring of X starting at index I. Use these two functions to produce a list of unique major site identifiers. (For this data, the list should contain only &quot;DR&quot; and &quot;MSK&quot;).</p>
</div>


<div class="keypoints">
<h4 id="key-points">Key Points</h4>
<ul>
<li>SQL can perform calculations using the values in a record as part of a query.</li>
</ul>
</div>
