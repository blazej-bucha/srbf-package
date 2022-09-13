subroutine DH_grid(nmax, lati, loni)
!
! =========================================================================
!
! DESCRIPTION: This function generates spherical latitudes and longitudes
!              of the Driscoll--Healy grid points on the unit sphere
!              (in radians) based on the input "nmax" value.
!
!
! INPUT: "nmax" -- Maximum harmonic degree; a non-negative integer
!
!
! OUTPUTS: "lati" -- Spherical latitudes (in radians) of the Driscoll--Healy
!                    grid points at a single meridian; a vector of dimensions
!                    (2 * (nmax + 1) + 1, 1)
!
!          "loni" -- Spherical longitudes (in radians) of the Driscoll--Healy
!                    grid points at a single parallel, a vector of dimensions
!                    (2 * (nmax + 1), 1);
!
!
! REFERENCES: Driscoll, J. R., Healy, D.M., 1994. Computing Fourier transforms
!                and convolutions on the 2-sphere. Advances in Applied
!                Mathematics 15, 202-250, doi:
!                http://doi.org/10.1006/aama.1994.1008
!
!             Schmidt, M., Fengler, M., Mayer-Gurr, T., Eicker, A., Kusche,
!                J., 2007. Regional gravity modelling in terms of spherical
!                base functions. Journal of Geodesy 81, 17-38,
!                doi: https://doi.org/10.1007/s00190-006-0101-5
!
!
! Contact: blazej.bucha@stuba.sk
!
!
! Code history: Version 1.0 (Feb 26, 2020)
!
!                           -- The first published version of the code
!
! =========================================================================


    ! Import modules
    ! ===============================================
    use vartypes  ! Module defining floating point double precision numbers
    use constants ! Module defining constants (pi, etc.)
    ! ===============================================


    ! Declaration of variables
    ! ===============================================
    implicit none


    integer,  intent(in)  :: nmax                     ! Maximum harmonic degree
    real(dp), intent(out) :: lati(2 * (nmax + 1) + 1) ! Latitudes of the
                                                      ! Driscoll--Healy grid
    real(dp), intent(out) :: loni(2 * (nmax + 1))     ! Longitudes of the
                                                      ! Driscoll--Healy grid
    integer :: i, l ! Useful substitutions
    ! ===============================================


    ! The main computational part of the routine
    ! ===============================================
    l = nmax + 1

    ! Compute the latitudes
    do i = 0, (2 * l)
        lati(i + 1) = -pi / 2.0_dp + pi * dble(i) / (2.0_dp * dble(l))
    end do

    ! Compute the longitudes
    do i = 0, (2 * l - 1)
        loni(i + 1) = pi * dble(i) / dble(l)
    end do
    ! ===============================================

end subroutine DH_grid
