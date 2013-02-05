# pass_array_args.py
import numpy as np
import _array_args

print _array_args.array_args.__doc__

# int_arr is a 10 X 10 array filled with consecutive integers.
# It is in 'fortran' order.
int_arr = np.asfortranarray(np.arange(100, dtype = 'i').reshape(10,10))

# cplx_arr is a 10 X 10 complex array filled with zeros.
# It is in 'fortran' order.
cplx_arr = np.asfortranarray(np.zeros((10,10), dtype = 'F'))

# We invoke the wrapped fortran subroutine.
real_arr = _array_args.array_args(int_arr, cplx_arr)

# Here are the results.
print "int_arr  = %s" %  int_arr
print "real_arr = %s" % real_arr
print "cplx_arr = %s" % cplx_arr
