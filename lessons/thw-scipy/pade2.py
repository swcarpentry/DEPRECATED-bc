#The Hacker Within: Python Boot Camp 2010 - Session 07 - Using SciPy.
#Presented by Anthony Scopatz.
#
#SciPy Pade, glide before you fly!

#As you have seen, SciPy has some really neat functionality that comes stock.
#Oddly, some of the best stuff is in the 'miscelaneous' module.
import scipy.misc 
from pylab import *

#So our exponential pade approimation didn't give us great gains, 
#But let's try approximating a rougher function.
def f(x):
	return (7.0 + (1+x)**(4.0/3.0))**(1.0/3.0)

#Through someone else's labors we know the expansion to be... 
f_exp = [2.0, 1.0/9.0, 1.0/81.0, -49.0/8748.0, 175.0/78732.0]

#The Pade coefficients are given simply by, 
p, q = scipy.misc.pade(f_exp, (5-1)/2)
#p and q are of numpy's polynomial class
#So the Pade approximation is given by 
def PadeAppx(x):
	return p(x) / q(x)

#Let's test it...
x = arange(0.0, 10.01, 0.01)

f_exp.reverse()
f_poly = poly1d(f_exp)

plot(x, PadeAppx(x), 'k--', label="Pade Approximation")
plot(x, f(x), 'k-', label=r'$f(x)$')
plot(x, f_poly(x), 'r-', label="Power Series")

xlabel(r'$x$')
ylabel("Polynomial Function")

legend(loc=0)

show()

#Check out http://docs.scipy.org/doc/scipy/reference/misc.html for a complete listing.

