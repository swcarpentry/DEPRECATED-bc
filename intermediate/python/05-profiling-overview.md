---
layout: lesson
root: ../..
---


# Introduction to Profiling in IPython Notebook


<pre class="in"><code>%matplotlib inline</code></pre>


## What is Profiling?

- "dynamic programme analysis to determine resource requirements" -- Wikipedia
- code benchmarking get metrics like runtime and memory usage



## Why Profile?


    "Premature optimization is the root of all evil." -- Knuth


- code should first be written for *clarity*
    - program correctness is more important than program speed
    - programmer time is more expensive than machine time


- sometimes speed is important, BUT
    - the slowdown might not be where you think it is
    - "improvements" might or might not be helpful
    - ugly code requires quantitative empirical justification


## Optimization Workflow (from SciPy Lecture notes 2.4.2)

1. make it work -- write code that can be understood
2. make it work reliably -- add test cases
3. use profiling to identify code bottlenecks and focus on these
    - look for the 10% effort that will give 90% benefit
    


## iPython Magics

- special commands that can be run using "%magic" in an ipython window
- built-in profiling magics:
    - time
    - timeit
    - prun



## time magic

- simplest approach, lowest overhead
- only runs code *once*
- outputs CPU times: user process, system process, total, and "Wall time" (Wall clock / stopwatch time)


<pre class="in"><code>def fib(n):
    return n if n &lt; 2 else fib(n - 1) + fib(n - 2)

%time fib(40)</code></pre>

<div class="out"><pre class='out'><code>CPU times: user 31 s, sys: 86.6 ms, total: 31.1 s
Wall time: 31 s
</code></pre><pre class='out'><code>102334155</code></pre></div>


## timeit magic

- runs code repeatedly (10,000 times by default) and measures total time required
- repeats this 3 times and reports on best
- best for timing small/short segments


<pre class="in"><code>def fib(n):
    return n if n &lt; 2 else fib(n - 1) + fib(n - 2)

# same as &#34;import timeit; timeit.timeit(&#39;fib(10)&#39;)&#34;
%timeit fib(10)</code></pre>

<div class="out"><pre class='out'><code>100000 loops, best of 3: 16.3 µs per loop
</code></pre></div>


## prun magic

- uses cProfile module to collect statistics
- tracks how many times functions are called, how long each call takes
    - helps determine which functions are running slowly
    - helps determine which functions are important to optimize


<pre class="in"><code>from time import sleep

def foo(): 
    sleep(1)

def bar(): 
    sleep(2)

def baz(): 
    foo()
    bar()

#same as &#34;import cProfile; cProfile.run(&#39;baz()&#39;)
stats = %prun -D baz.prof -qr baz()

stats.sort_stats(-1).print_stats()</code></pre>

<div class="out"><pre class='out'><code> 
*** Profile stats marshalled to file u&#39;baz.prof&#39;. 
         7 function calls in 3.003 seconds

   Ordered by: standard name

   ncalls  tottime  percall  cumtime  percall filename:lineno(function)
        1    0.000    0.000    1.001    1.001 &lt;ipython-input-4-11141aa6748f&gt;:3(foo)
        1    0.000    0.000    2.002    2.002 &lt;ipython-input-4-11141aa6748f&gt;:6(bar)
        1    0.000    0.000    3.003    3.003 &lt;ipython-input-4-11141aa6748f&gt;:9(baz)
        1    0.000    0.000    3.003    3.003 &lt;string&gt;:1(&lt;module&gt;)
        1    0.000    0.000    0.000    0.000 {method &#39;disable&#39; of &#39;_lsprof.Profiler&#39; objects}
        2    3.003    1.502    3.003    1.502 {time.sleep}


</code></pre><pre class='out'><code>&lt;pstats.Stats instance at 0x7feee2d29c20&gt;</code></pre></div>


## Understanding cProfile Output

*ncalls*
:    the number of times the function is called

*tottime*
:    the total time spent executing the function across multiple calls, *excluding* calls to other functions

*percall*
:    tottime / ncalls

*cumtime*
:    total time spent in the function, *including* calls to subfunctions; accurate even for recursive functions

*percall* (2nd column)
:    cumtime / "primitive calls" (calls not caused by recursion)



## RunSnakeRun

- standard tool for visualisation
- easy to use: "runsnakerun baz.prof"



## Additional Extensions

These can all help out with profiling, but currently require a bit of time and effort to get installed correctly.

- lprun
- mprun
- memit
- snakeviz


## References

- http://scipy-lectures.github.io/ (particularly section 2.4.2)
- http://pynash.org/2013/03/06/timing-and-profiling.html
- http://www.huyng.com/posts/python-performance-analysis/
- http://www.vrplumber.com/programming/runsnakerun/
