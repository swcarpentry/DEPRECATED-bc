import numpy as np
import pylab
import pdb


#read in the image
img = np.genfromtxt('spec_example.dat')
ax_img = pylab.axes([0.1, 0.1, 0.65, 0.8]) #[left, bottom, width, height]
ax_plot = pylab.axes([0.77, 0.1, 0.13, 0.8])




#Display the image
ax_img.imshow(img, origin = 'lower', interpolation = 'nearest')

#Collapse the spectrum along x axis
img_collapse = np.sum(img, axis = 1)
#create and array to plot against
y = np.arange(img_collapse.shape[0])

#Plot to new axis
ax_plot.plot(img_collapse, y, 'k', lw = 2)
ax_plot.set_ylim(ax_img.get_ylim())
