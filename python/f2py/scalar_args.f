      subroutine scalar_args(int_in, real_in, int_inout, real_inout,
     \ int_out, real_out)
C Here are the f2py-specific comments.
Cf2py intent(in) :: int_in, real_in }}
Cf2py intent(inout) :: int_inout, real_inout
Cf2py intent(out) :: int_out, real_out

      int_inout = int_in
      real_inout = real_in
      int_out = int_inout
      real_out = real_inout
        
      end subroutine scalar_args
