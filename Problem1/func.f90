module func
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
          implicit none
          private
          public :: set_mode_2p, set_mode_3s, run_func, run_func_derv
          integer :: mode = 0          ! 0 unset, 1 = 2p orbital, 2 = 3s orbital

contains

          !Switch the mode for the orbital at hand
          subroutine set_mode_2p()
                    mode = 1
          end subroutine set_mode_2p

          subroutine set_mode_3s()
                    mode = 2
          end subroutine set_mode_3s

          pure function run_func(muR) result(val)
                    ! f(muR) = e^{-muR} * P(muR) - 0.05   (Eq. 8 / Eq. 10 rearranged)
                    real(dp), intent(in) :: muR
                    real(dp) :: val, ex, poly

                    ex = exp(-muR)

                    select case (mode)
                    case (1)   ! 2p: P = 1 + x + x^2/2 + x^3/6 + x^4/24
                              poly = 1.0_dp + muR + muR**2/2.0_dp + muR**3/6.0_dp + muR**4/24.0_dp
                              val = ex*poly - 0.05_dp
                    case (2)   ! 3s: P = 1 + x + x^2/2 + x^4/4 - x^5/12 + x^6/72
                              poly = 1.0_dp + muR + muR**2/2.0_dp + muR**4/4.0_dp - muR**5/12.0_dp + muR**6/72.0_dp
                              val = ex*poly - 0.05_dp
                    case default
                              val = ieee_value(val, ieee_quiet_nan)
                    end select

          end function run_func

          pure function run_func_derv(muR) result(val)
                    ! f'(muR) = e^{-muR} * (P' - P)   (product rule)
                    real(dp), intent(in) :: muR
                    real(dp) :: val, ex, dpoly

                    ex = exp(-muR)

                    select case (mode)
                    case (1)   ! 2p: P' - P = -x^4/24
                              dpoly = -muR**4/24.0_dp
                              val = ex*dpoly
                    case (2)   ! 3s: P' - P = -(x^2/2 - x^3 + 2x^4/3 - x^5/6 + x^6/72)
                              dpoly = -(muR**2/2.0_dp - muR**3 + 2.0_dp*muR**4/3.0_dp - muR**5/6.0_dp + muR**6/72.0_dp)
                              val = ex*dpoly
                    case default
                              val = ieee_value(val, ieee_quiet_nan)
                    end select

          end function run_func_derv

end module func
