# -*- coding: utf-8 -*-
# <nbformat>3</nbformat>

# <markdowncell>

# # Exercise 1
# 1. Assign variables with the values 1, 5, and "ten" with the types int, float, and string respectively.
# <codecell>
  
# <markdowncell>
# 2. Confirm that the types are int, float, and string
# <codecell>
  
# <markdowncell>
# 3. Determine which for which pairs of the set of int, float, and string the + operation gives an error.  
# > int-int :
# > int-float : 
# > int-string :
# > float-float  :
# > float-string :
# > string-string : 
# Any surprises?
# <codecell>
 
# <markdowncell>
# 4. Determine which for which pairs of the set of int, float, and string) the * operation gives an error.  Any surprises?
# <codecell>
  
# <markdowncell>
# 5. Assign a string the value of "1, 5, and ten" from these three variables.  
# <codecell>
 
# <markdowncell>
# # Exercise 2
# Here you will use **math.log10()** and **math.floor()**, which require the line **import math** for you to access these funcitons.  
# 1. Determine the return type of log10()
# <codecell>
import math
  
# <markdowncell>
# 2. What is the value and type of log10(42) ?
# <codecell>
  
# <markdowncell>
# 3. What is the value and type of log10(-0.32) ?
# <codecell>
 
# <markdowncell>
# 4. What about 1.0 / 0 ?   
# <codecell>
   
# <markdowncell>
# 4. What is the return type of floor acting on an int?  Acting on a float?  **floor** is in the math namespace.  It will only work after **import math** and it is invoked as **math.floor()**
# <codecell>
   
# <markdowncell>
# # Exercise 3
#  len() is a builtin function to count the length of things.  For which of the basic datatypes so far does len() return a value?  Does it return the length or the length+1 ?
# <codecell>
  
# <markdowncell>
# # Example 1
# Python lists are agnostic to the type of the data contained within them.  You can generate arrays with different elements of different types:
# <codecell>
pythonlist =[2.43, 1, 0.92, 0, "0.38"]
print pythonlist
# <codecell>
print type(pythonlist[0]), type(pythonlist[1]), type(pythonlist[2]), type(pythonlist[3])
# <markdowncell>
# # Exercise 4
# numpy is an extremely useful library that has its own data types; you will frequently need to specify the types as float or possibly higher-precision float at the time you create them.  The numpy data types are only sometimes compatible with the native data types.  
# <codecell>
import numpy as np
numpyarray = np.array(pythonlist)
# <markdowncell>
# What is the default data type of a after the conversion?
# <codecell>
 
# <markdowncell>
# Ack! The results of type() do not betray the representation.  For that we need 
# <codecell>
print numpyarray.dtype
# <markdowncell>
# Which is not a numeric data type.  We can cause it to be stored as numpy floats if we specify float when we convert it to numpy:
# <codecell>
numpyfloatarray = np.array(python_array, dtype="float")

# <markdowncell>
# # Exercise 5
# 5A. Write an expression to determine the number of digits in a non-negative integer.  Hint:  maybe **len()** or **math.log()** might be useful here.
# <codecell>

# <markdowncell>
# 5B. Test your expression on 45, 2, and 2.0.  Does it work for all three?
# <codecell>
