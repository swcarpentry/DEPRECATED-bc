#The Hacker Within: Python Boot Camp 2010 - Session 07 - Using SciPy.
#Presented by Anthony Scopatz.
#
#SciPy Image Tricks, fly before you....You can do that?!

#For some reason that has yet to be explained to me, SciPy has the ability to treat 2D & 3D arrays
#as images.  You can even convert PIL images or read in external files as numpy arrays!
#From here, you can fool around with the raw image data at will.  Naturally, this functionality 
#is buried within the 'miscellaneous' module.
import scipy.misc

#First let's read in an image file.  For now, make it a JPEG.
img = scipy.misc.imread("image.jpg")
#Note that this really is an array!
print(str(img))

#We can now apply some basic filters...
img = scipy.misc.imfilter(img, 'blur')

#We can even rotate the image, counter-clockwise by degrees.
img = scipy.misc.imrotate(img, 45)

#And then, we can rewrite the array to an image file.
scipy.misc.imsave("image1.jpg", img)

#Because the array takes integer values from 0 - 255, we can easily define our own filters as well!
def InverseImage(imgarr):
	return 255 - imgarr

#Starting fresh we get... 
img = scipy.misc.imread("image.jpg")
img = scipy.misc.imrotate(img, 330)
img = InverseImage(img)
scipy.misc.imsave("image2.jpg", img)

#Check out http://docs.scipy.org/doc/scipy/reference/misc.html for a complete listing.

