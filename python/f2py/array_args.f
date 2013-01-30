      subroutine array_args(nx, ny, int_arr_in,
     \                      cmplx_arr_inout, 
     \                      real_arr_out)

          integer nx, ny
          integer int_arr_in(nx, ny)
          complex cmplx_arr_inout(nx, ny)
          real real_arr_out(nx, ny)

Cf2py intent(in) nx, ny
Cf2py intent(in) int_arr_in
Cf2py intent(inout) cmplx_arr_inout
Cf2py intent(out) real_arr_out

C ... body of subroutine unchanged ...

        integer i, j

        do j = 1, ny
            do i = 1, nx
                cmplx_arr_inout(i,j) = cmplx(int_arr_in(i,j),
   \                   int_arr_in(i,j))
                real_arr_out(i,j) = real(int_arr_in(i,j))
            enddo
        enddo

      end subroutine array_args
