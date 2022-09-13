subroutine sph2cart(lon, lat, n, x, y, z)
!
! ==========================================================================
!
! DESCRIPTION: This function transforms spherical coordinates into cartesian
!              coordinates, assuming that the points reside on the unit sphere.
!
!
! INPUT: "lon" -- vector of spherical longitudes in radians
!        "lat" -- vector of spherical latitudes in radians
!        "n"   -- number of input points
!
!
! OUTPUT: "x", "y", "z" -- three vectors of cartesian coordinates
!
!
! Contact: blazej.bucha@stuba.sk
!
!
! Code history: Version 1.0 (Feb 26, 2020)
!
!                           -- The first published version of the code
!
! ==========================================================================


    ! Import modules
    ! ===============================================
    use vartypes  ! Module defining floating point double precision numbers
    use constants ! Module defining constants (pi, etc.)
    ! ===============================================


    ! Declaration of variables
    ! ===============================================
    implicit none


    integer,  intent(in)  :: n                ! Number of points
    real(dp), intent(in)  :: lon(n), lat(n)   ! Input longitudes and latitudes
    real(dp), intent(out) :: x(n), y(n), z(n) ! Output cartesian coordinates
    integer               :: i

    ! ===============================================


    ! The main computational part of the routine
    ! ===============================================
    do i = 1, n

        x(i) = cos(lat(i)) * cos(lon(i))
        y(i) = cos(lat(i)) * sin(lon(i))
        z(i) = sin(lat(i))

    end do
    ! ===============================================

end subroutine sph2cart
