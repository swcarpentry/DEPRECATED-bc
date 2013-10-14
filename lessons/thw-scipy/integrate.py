#The Hacker Within: Python Boot Camp 2010 - Session 07 - Using SciPy.
#Presented by Anthony Scopatz.
#
#SciPy Integration, run before you glide. 

#Tools used to calculate numerical, definite integrals may be found in the 'integrate' module.
import scipy.integrate
#For kicks, let's also grab
import scipy.special
import numpy

#There are two basic ways you can integrate in SciPy:
#     1. Integrate a function, or
#     2. Integrate piecewise data.

#First Let's deal with integration of functions.
#Recall that in Python, functions are also objects.  
#Therefore you can pass functions as arguments to other functions!
#Just make sure that the function that you want to integrate returns a float, 
#or, at the very least, an object that has a __float__() method.

#The simplest way to compute a functions definite integral is via the quad(...) function.
def CrazyFunc(x):
	return (scipy.special.i1(x) - 1)**3

print("Try integrating CrazyFunc on the range [-5, 10]...")

val, err = scipy.integrate.quad(CrazyFunc, -5, 10)

print("A Crazy Function integrates to %.8E"%val)  
print("And with insanely low error of %.8E"%err)  
print("")

#You can also use scipy.integrate.Inf for infinity in the limits of integration
print("Now try integrating e^x on [-inf, 0]")
print("(val, err) = " + str( scipy.integrate.quad(scipy.exp, -scipy.integrate.Inf, 0.0) ))
print("")

#2D integrations follows similarly, 
def dA_Sphere(phi, theta):
	return  scipy.sin(phi)

print("Integrate the surface area of the unit sphere...")
val, err = scipy.integrate.dblquad(dA_Sphere, 0.0, 2.0*scipy.pi, lambda theta: 0.0,  lambda theta: scipy.pi )
print("val = %.8F"%val)
print("err = %.8E"%err)
print("")

def dV_Sphere(phi, theta, r):
	return r * r * dA_Sphere(phi, theta)

print("Integrate the volume of a sphere with r=3.5...")
val, err = scipy.integrate.tplquad(dV_Sphere, 0.0, 3.5, lambda r: 0.0, lambda r: 2.0*scipy.pi, lambda x, y: 0.0, lambda x, y: scipy.pi)
print("val = %.8F"%val)
print("err = %.8E"%err)
print("")

#Now, only very rarely will scientists (and even more rarely engineers) will truely 'know' 
#the function that they wish to integrate.  Much more often we'll have piecewise data 
#that we wish numerically integrate (ie sum an array y(x), biased by array x).  
#This can be done in SciPy through the trapz function.

y = range(0, 11)
print("Trapazoidally integrate y = x on [0,10]...")
val = scipy.integrate.trapz(y)
print("val = %F"%val)
print("")

#You can also define a domain to integrate over.
x = numpy.arange(0.0, 20.5, 0.5)
y = x * x
print("Trapazoidally integrate y = x^2 on [0,20] with half steps...")
val = scipy.integrate.trapz(y, x)
print("val = %F"%val)
print("")

print("Trapazoidally integrate y = x^2 with dx=0.5...")
val = scipy.integrate.trapz(y, dx=0.5)
print("val = %F"%val)
print("")

def dDecay(y, t, lam):
	return -lam*y

#Of course, sometimes we have simple ODEs that we want to integrate over time for...
#These are generally of the form:
#     dy / dt = f(y, t)
#For example take the decay equation...
#     f(y, t) = - lambda * y
#We can integrate this using SciPy's 'odeint'  This is of the form:
#     odeint( f, y0, [t0, t1, ...])
#Let's try it... 
vals = scipy.integrate.odeint( lambda y, t: dDecay(y, t, 0.2), 1.0, [0.0, 10.0] ) 
print("If you start with a mass of y(0) = %F"%vals[0][0])
print("you'll only have y(t=10) = %F left."%vals[1][0])

#Check out http://docs.scipy.org/doc/scipy/reference/integrate.html for a complete listing.
