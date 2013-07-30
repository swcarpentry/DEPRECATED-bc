#The Hacker Within: Python Boot Camp 2010 - Session 07 - Using SciPy.
#Presented by Anthony Scopatz.
#
#SciPy special functions, walk before you run!

#Code that numerically approximates common (and some not-so-common) special functions can be found in 'scipy.special'
from scipy.special import *

#Here you can find things like error functions, gamma functions, Legendre polynomials, etc.
#But as a example let's focus on my favorites: Bessel functions.
#Time for some graphs...
from pylab import *

x = arange(0.0, 10.1, 0.1)

for n in range(4):
	j = jn(n, x)
	plot(x, j, 'k-')
	text(x[10*(n+1)+1], j[10*(n+1)], r'$J_%r$'%n)

for n in range(3):
	y = yn(n, x)
	plot(x, y, 'k--')
	text(x[10*(n)+6], y[10*(n)+5], r'$Y_%r$'%n)

axis([0, 10, -2, 1.25])
xlabel(r'$x$')
ylabel("Bessel Functions")

show()

#Check out http://docs.scipy.org/doc/scipy/reference/special.html for a complete listing.

#Note that the figure that was created here is a reproduction of 
#Figure 6.5.1 in 'Numerical Recipes' by W. H. Press, et al.
