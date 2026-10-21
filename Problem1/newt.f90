!Newtons Method for Calcuating Roots - Zackary Childs Sept 15, 2025

module newt
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
          use func, only: run_func, run_func_derv
          implicit none
          private
          public :: newton

contains

          pure subroutine newton(a_in, tol_p, root, froot, iters, status, maxit)

                    real(dp), intent(in) :: a_in, tol_p
                    real(dp), intent(out) :: root, froot
                    integer, intent(out) :: iters, status
                    !For status 0 = Ok, 1 = Safe Stop, 2 = Max Iterations hit
                    integer, intent(in), optional :: maxit
                    real(dp) :: TOL, p, p0, fp0, fp_p0
                    integer :: kmax, k
                    real(dp), parameter :: safeden = 100.0_dp*epsilon(1.0_dp)

                    TOL = max(epsilon(1.0_dp), merge(tol_p, 1.0e-12_dp, tol_p > 0.0_dp))

                    !Check for maxit otherwise create one based on the b and a values
                    if (present(maxit)) then
                              kmax = maxit
                    else
                              kmax = 2000000
                    end if

                    p0 = a_in

                    !Actual Newtons Method
                    do k = 1, kmax
                              fp0 = run_func(p0)
                              fp_p0 = run_func_derv(p0)

                              if (abs(fp_p0) <= safeden*(1.0_dp + abs(p0))) then
                                        root = p0
                                        froot = fp0
                                        iters = k - 1
                                        status = 1
                                        return
                              end if

                              p = p0 - (fp0/fp_p0)

                              !Catch overflow / NaN instead of iterating on garbage
                              if (.not. ieee_is_finite(p)) then
                                        root = p0
                                        froot = fp0
                                        iters = k
                                        status = 3
                                        return
                              end if

                              if (abs(p - p0) < TOL) then
                                        root = p
                                        froot = run_func(p)
                                        iters = k
                                        status = 0
                                        return
                              else
                                        p0 = p
                              end if

                    end do
                    !What happens if method does not find the root after the maximum number of iterations
                    root = p0
                    froot = run_func(root)
                    iters = kmax
                    status = 2

          end subroutine newton

end module newt

