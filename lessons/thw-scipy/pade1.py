#The Hacker Within: Python Boot Camp 2010 - Session 07 - Using SciPy.
#Presented by Anthony Scopatz.
#
#SciPy Pade, glide before you fly!

#As you have seen, SciPy has some really neat functionality that comes stock.
#Oddly, some of the best stuff is in the 'miscelaneous' module.
import scipy.misc 

#Most people are familar with the polynomial expansions of a function:
#     f(x) = a + bx + cx^2 + ...
#Or a Taylor expansion:
#     f(x) = sum( d^n f(a) / dx^n (x-a)^n /n! )
#However, there exists the lesser known, more exact Pade approximation.
#This basically splits up a function into a numerator and a denominator.
#     f(x) = p(x) / q(x)
#Then, you can approximate p(x) and q(x) using a power series.  
#A more complete treatment is available in Section 5.12 in 'Numerical Recipes' by W. H. Press, et al.


#The stregnth of this method is demonstated though figures...
import numpy as np 
from matplotlib import pyplot as plt 

#Let's expand e^x to fith order and record the coefficents 
e_exp = [1.0, 1.0, 1.0/2.0, 1.0/6.0, 1.0/24.0, 1.0/120.0]

#The Pade coefficients are given simply by, 
p, q = scipy.misc.pade(e_exp, 2)
#p and q are of numpy's polynomial class
#So the Pade approximation is given by 
def PadeAppx(x):
	return p(x) / q(x)

#Let's test it...
x = np.arange(0.0, 3.1, 0.1)

e_exp.reverse()
e_poly = np.poly1d(e_exp)

plt.plot(x, PadeAppx(x), 'k--', label="Pade Approximation")
plt.plot(x, scipy.e**x, 'k-', label=r'$e^x$')
plt.plot(x, e_poly(x), 'r-', label="Power Series")

#axis([0, 10, -2, 1.25])
plt.xlabel(r'$x$')
plt.ylabel("Exponential Functions")

plt.legend(loc=0)

plt.show()

#Check out http://docs.scipy.org/doc/scipy/reference/misc.html for a complete listing.
