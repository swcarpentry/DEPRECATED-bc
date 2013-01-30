# pass_args.py
import numpy as np
import _scalar_args

print _scalar_args.scalar_args.__doc__

# these are simple python scalars.
int_in = 1.0
real_in = 10.0

# since these are intent(inout) variables, these must be arrays
int_inout = np.zeros((1,), dtype = np.int32)
real_inout = np.zeros((1,), dtype = np.float32)

# all intent(out) variables are returned in a tuple, so they aren't passed as
# arguments.

int_out, real_out = _scalar_args.scalar_args(int_in, real_in, int_inout, real_inout)

for name in ('int_inout', 'real_inout', 'int_out', 'real_out'):
    print '%s == %s' % (name, locals()[name])
