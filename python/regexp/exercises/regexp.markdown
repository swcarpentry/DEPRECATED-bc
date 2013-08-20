Please commit all the files to a folder named `regex` in your repository. 

## Warm up

1. How do we import the regular expressions library in Python?
   ```python    
   import re
   ```
2. We want to match `xx` or `xy` in a string with variable name `data`.  Why
will the following code not work properly?
    
   ```python    
    re.search(data, "xx|y")
   ```
   **Answer:**
   > 1. The order of the arguments is incorrect.  The first argument should be
   > the pattern, and the second should be the data.
   > 
   > 2. We have a problem with precedence: `xx|y` matches `xx` or `y`.  We need
   > to change it to `x(x|y)`.  (The pattern `xx|xy` would also work.)

3. Suppose we have a file called _people.txt_ that contains the following data:
    
    ```
    Rose	416-333-4444	rose@someplace.com
    Martha	905-888-1234	martha@hotmail.com
    Donna	647-222-9876	donna@rogers.ca
    Amy	905-777-2222	amy@gmail.com
    ```

   * Using regular expressions, write code that will read from the file, and
     print out the phone number for each individual.
     
     **Answer:**
     > ```python
     > import re
     > people = open('people.txt', 'r')
     > for p in people :
     >     match = re.search('...-...-....', p)
     >     print match.group(0)
     > ```

   * Write code that prints just the area code (the first three digits of a
     phone number) for each individual.
     
     **Answer:**
     > ```python
     > import re
     > people = open('people.txt', 'r')
     > for p in people:
     >     match = re.search('(...)-...-....', p)
     >     print match.group(1)
     > ```

## Operators

1. Give a pattern to match the following:
   * Words that consist entirely of uppercase letters, and contain 5 to 8
     letters.
            
            [A-Z]{5,8}
   * Words that begin with an uppercase letter, followed by any number
     (including zero!) of lowercase letters.

            [A-Z][a-z]*
   * One or more vowels, followed by an optional space character, followed by
     one or more numbers.

            [aeiou]+ ?[1-9]+

2. Write a function called time_match that takes in one parameter (a string).
Use regular expressions to match the string with the following two time formats:
   * HH:MM:SS in 24-hour format. Here, we require two digits for the hour,
     minutes, and seconds.
   * HH:MM PM/AM  in 12-hour format.  Here, we require one or two digits for
     the hour, two digits for the minutes, and the seconds should not be
     included.

   So, for example, "15:41:38" and "03:29:10" are valid matches for the first
   format, and "2:30 PM" and "10:54 AM" are valid matches for second format. 

   If the string matches one of the formats, the function should return the
   matching hours and minutes as "HH:MM".  For example, if the string contains
   "15:42:38", we should return "15:42".  If the string contains "2:30 PM", we
   should return "2:30".  

   If you are feeling ambitious, for the strings that match the second format,
   convert the result to 24-hour time before returning it.  (So for "2:30 PM" we
   would return "14:30".) 
   
   ```python
    def time_match(record):
        # HH:MM:SS
        m = re.search('([0-9]{2}):([0-9]{2}):[0-9]{2}', record)
        if m:
            return m.group(1) + ':' + m.group(2)    
    
        # (H)H:MM AM/PM
        m = re.search('([0-9]{1,2}):([0-9]{2}) (AM|PM)', record)
        if m:
            if (m.group(3) == 'AM'):
                return m.group(1) + ":" + m.group(2)
            else:
                # convert to 24-hour time format
                return str(int(m.group(1)) + 12) + ":" + m.group(2)
        return None
    ```


## Exercises: Cities of Canada

You are studying the propagation of an epidemic disease. Using the mathematical
models developed by your team, you wrote a program to simulate the spread of the
disease across several cities. For each city, your simulation requires the
population density and the location of the city (latitude and longitude).
Because you don't have this data at hand, you asked a colleage to send it to
you, but he hasn't responded yet. You would like to start testing your program
while you are waiting, so you head to Wikipedia and realize that the information
listed in http://en.wikipedia.org/wiki/List_of_cities_in_Canada seems to be good
enough for testing your algorithm. But first, you need to parse it...

### Exercise 1: inspecting the data

Open the wikipedia page for [List of cities in Canada - wikipedia
page](http://en.wikipedia.org/wiki/List_of_cities_in_Canada) and the [source
text for the wikipedia
page](http://svn.software-carpentry.org/swc/4.0/topics/regexp/exercises/List_of_cities_in_Canada.txt).
Some browsers may display the source text incorrectly. You may want to download
the text file and open it with a text editor instead of reading it in the
browser.

1. Can you see a pattern you could use to detect where in the source the section
for "Alberta" begins? Write a regular expression that matches this pattern.
Could you use or extend that pattern to find the sections for the other
provinces?. Save your answer in a file called `cities_patterns.txt` and commit
it to the repository.

2. The pattern you wrote, will it match only the section headers corresponding
to provinces, or will it match all section headers? Will it match all the
provinces? Write a program to test this out. Save the program a parse_cities.py.
You can use [this
module](http://svn.software-carpentry.org/swc/4.0/topics/regexp/exercises/wikipedia.py)
to download the source text from your program.


### Exercise 2: the first parser

Focus now on Alberta. The information for each city seems to be stored in one
line each and follow this format:

    || [[Airdrie, Alberta|Airdrie]] || 33.10|| {{nts|28927}} || {{nts|874}} || {{Coord|51.29173|114.0143|region:CA-AB_type:city_scale:50000|name=Airdrie}} || 

1. Write a regular expression to extract the relevant information (higlighted)
from this line. Test the regular expression in the python shell and append it to
the file `cities_patterns.txt`.

2. Write a program that prints the information for the cities in Alberta in a
more readable format. There are 17 cities in the wikipedia article - make sure
that your output includes the 17 cities. Note that the format for the city name
is not always the same, for instance, Calgary doesn't have a comma (`,`) or a
pipe (`|`).  Make sure to only print the cities in Alberta. Save the program in
the file `alberta_cities.py`. Hint: start parsing with the Alberta regular
expression after you find the Alberta section using the expression from the
previous exercise, and stop as soon as you find the next section.


### Exercise 3: other provinces

1. Will this regular expression produce the correct results for the cities in
Manitoba? Why? Append your answer to the file `cities_patterns.txt`

2. Write a function `cities()` to parse this page. The function should return a
list of dictionaries, containing the province, city name, population, area,
longitude and latitude of every city. For the purposes of this exercise, focus
only on the provinces of Alberta and Manitoba.  Save this function in a file
called `canadian_cities.py`. Include a `__main__` section to display this data.

**Hint:** Create a dictionary with the provinces as keys, and the regular
expressions and order of parameters as the values, like it was done in the
episode http://software-carpentry.org/4_0/regexp/patterns/.  When you find the
section of a province that you have in the dictionary, start using the
corresponding regular expression until you find the next section.
