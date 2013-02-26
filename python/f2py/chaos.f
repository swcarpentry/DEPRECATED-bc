    subroutine iterate_limit(func, x0, num_iters, results, n)
        external func
        double precision func
        double precision x0
        integer num_iters, n
        double precision results(n)

        integer i

        do i = 1, num_iters
            x0 = func(x0)
        enddo

        do i = 1, n
            results(i) = x0
            x0 = func(x0)
        enddo

    end subroutine iterate_limit
