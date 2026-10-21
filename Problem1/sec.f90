!Secants Method for Calculating Root - Zackary Childs Sept 15, 2025

module sec
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
          use func, only: run_func
          implicit none
          private
          public :: secant

contains

          pure subroutine secant(a_in, b_in, tol_p, root, froot, iters, status, maxit)

                    real(dp), intent(in) :: a_in, b_in, tol_p
                    real(dp), intent(out) :: root, froot
                    integer, intent(out) :: iters, status
                    !For status 0 = Ok, 1 = Max Iterations hit
                    integer, intent(in), optional :: maxit
                    real(dp) :: TOL, p1, p0, fp1, fp0, p
                    integer :: kmax, k
                    real(dp), parameter :: safeden = 100.0_dp*epsilon(1.0_dp)

                    TOL = max(epsilon(1.0_dp), merge(tol_p, 1.0e-12_dp, tol_p > 0.0_dp))

                    !Check for maxit otherwise create one based on the b and a values
                    if (present(maxit)) then
                              kmax = maxit
                    else
                              kmax = 2000000
                    end if

                    !Ensure ordered interval
                    if (a_in <= b_in) then
                              p0 = a_in
                              p1 = b_in
                    else
                              p0 = b_in
                              p1 = a_in
                    end if

                    fp0 = run_func(p0)
                    fp1 = run_func(p1)

                    !Actual Secant Method
                    do k = 2, kmax

                              !Flat secant line: the step would divide by ~0
                              if (abs(fp1 - fp0) <= safeden*max(abs(fp1), abs(fp0))) then
                                        root = p1
                                        froot = fp1
                                        iters = k
                                        status = 1
                                        return
                              end if

                              p = p1 - ((fp1*(p1 - p0))/(fp1 - fp0))

                              !Catch overflow / NaN instead of iterating on garbage
                              if (.not. ieee_is_finite(p)) then
                                        root = p1
                                        froot = fp1
                                        iters = k
                                        status = 3
                                        return
                              end if

                              if (abs(p - p1) < TOL) then
                                        root = p
                                        froot = run_func(root)
                                        iters = k
                                        status = 0
                                        return
                              end if

                              p0 = p1
                              fp0 = fp1
                              p1 = p
                              fp1 = run_func(p)

                    end do

                    !What happens if method does not find the root after the maximum number of iterations

                    p = p1 - ((fp1*(p1 - p0))/(fp1 - fp0))
                    root = p
                    froot = run_func(root)
                    iters = kmax
                    status = 2

          end subroutine secant

end module sec

