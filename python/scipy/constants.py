#The Hacker Within: Python Boot Camp 2010 - Session 07 - Using SciPy.
#Presented by Anthony Scopatz.
#
#SciPy constants, Crawl before you walk!

#A plethora of important fundamental constants can be found in
import scipy.constants
#NOTE: this module is not automatically included when you "import scipy"

#Some very basic pieces of information are given as module attributes
print("SciPy thinks that pi = %.16f"%scipy.constants.pi)
import math
print("While math thinks that pi = %.16f"%math.pi)
print("SciPy also thinks that the speed of light is c = %.1F"%scipy.constants.c)
print("")

#But the real value of SciPy constants is its enormous physical constant database
print("SciPy physical constants are of the form:")
print("      scipy.constants.physical_constants[name] = (value, units, uncertainty)")
print("")

print("For example the mass of an alpha particle is %s"%str(scipy.constants.physical_constants["alpha particle mass"]))
print("But buyer beware! Let's look at the speed of light again.")
print("c = %s"%str(scipy.constants.physical_constants["speed of light in vacuum"]))
print("The uncertainty in c should not be zero!")
print("")

print("Check http://docs.scipy.org/doc/scipy/reference/constants.html for a complete listing.")

