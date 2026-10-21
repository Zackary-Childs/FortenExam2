program Radius95
          use, intrinsic :: iso_fortran_env, only: dp => real64
          use newt
          use sec
          use func
          implicit none

          real(dp), parameter :: tol = 1.0e-10_dp
          real(dp), parameter :: IE_H = 13.6_dp        ! eV, hydrogen 1s (Hayward Eq. 7.37)
          real(dp), parameter :: a0 = 0.529177_dp    ! Bohr radius in Angstrom

          !Atom data: name, principal quantum number of outer electron, NIST first ionization energy (eV)
          character(len=2), parameter :: atom(3) = ['O ', 'Ne', 'Mg']
          integer, parameter :: n(3) = [2, 2, 3]
          real(dp), parameter :: IE(3) = [13.618_dp, 21.565_dp, 7.646_dp]

          real(dp) :: muR2p, muR3s, root, froot
          real(dp) :: Zeff(3), mu(3), R(3)
          integer :: iters, status, i, j, tmp
          integer :: order(3) = [1, 2, 3]

          !2p orbital (O, Ne): solve Eq. 8 for muR
          call set_mode_2p()

          call newton(a_in=9.0_dp, tol_p=tol, root=muR2p, froot=froot, iters=iters, status=status)
          call report('2p Newton', muR2p, froot, iters, status)

          call secant(a_in=9.0_dp, b_in=10.0_dp, tol_p=tol, root=root, froot=froot, iters=iters, status=status)
          call report('2p Secant', root, froot, iters, status)

          !3s orbital (Mg): solve Eq. 10 for muR
          call set_mode_3s()

          call newton(a_in=14.0_dp, tol_p=tol, root=muR3s, froot=froot, iters=iters, status=status)
          call report('3s Newton', muR3s, froot, iters, status)

          call secant(a_in=14.0_dp, b_in=15.0_dp, tol_p=tol, root=root, froot=froot, iters=iters, status=status)
          call report('3s Secant', root, froot, iters, status)

          !Convert muR to R for each atom: Zeff = n*sqrt(IE/13.6), mu = 2*Zeff/n, R = muR/mu
          print *, ''
          print '(a)', '  Atom  Orbital   Zeff     R95 (a0)   R95 (A)'
          do i = 1, 3
                    Zeff(i) = n(i)*sqrt(IE(i)/IE_H)
                    mu(i) = 2.0_dp*Zeff(i)/n(i)
                    R(i) = merge(muR2p, muR3s, n(i) == 2)/mu(i)
                    print '(2x,a2,4x,a2,6x,f6.4,3x,f8.4,3x,f7.4)', atom(i), merge('2p', '3s', n(i) == 2), Zeff(i), R(i), R(i)*a0
          end do

          !Rank largest to smallest
          do i = 1, 2
                    do j = i + 1, 3
                              if (R(order(j)) > R(order(i))) then
                                        tmp = order(i); order(i) = order(j); order(j) = tmp
                              end if
                    end do
          end do

          print *, ''
          print '(a,a,a,a,a,a)', '  Ranking: ', trim(atom(order(1))), ' > ', trim(atom(order(2))), ' > ', trim(atom(order(3)))

contains

          !Print one root-finder result and stop if it did not converge
          subroutine report(label, x, fx, it, st)
                    character(len=*), intent(in) :: label
                    real(dp), intent(in) :: x, fx
                    integer, intent(in) :: it, st

                    print '(2x,a,a,f12.8,a,i0,a,es10.3,a,i0)', label, ':  muR = ', x, '  iters = ', it, &
                              '  |f| = ', abs(fx), '  status = ', st

                    if (st /= 0) error stop 'Root finder did not converge -- check the starting guess'

          end subroutine report

end program Radius95
