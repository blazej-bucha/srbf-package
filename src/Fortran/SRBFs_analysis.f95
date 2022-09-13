subroutine SRBFs_analysis(lativ, loniv, t, radius, nmax, lati, loni, a)
!
! =========================================================================
!
! DESCRIPTION: This function performs surface SRBF analysis up to degree "nmax"
!              of a function sampled at the Driscoll--Healy grid residing
!              on the reference sphere with the radius "radius".
!
!              The loop over latitude parallels is parallelized using OpenMP.
!
!
! INPUTS: "lativ"  -- Vector of spherical latitudes (in radians) over a single
!                     meridian of the Driscoll--Healy grid (the first output
!                     from the "DH_grid" function); dimensions of the vector
!                     (2 * (nmax + 1) + 1, 1)
!
!         "loniv"  -- Vector of spherical longitudes (in radians) over a single
!                     parallel of the Driscoll--Healy grid (the second output
!                     from the "DH_grid" function); dimensions of the vector
!                     (2 * (nmax + 1), 1)
!
!         "t"      -- Input data that are analysed; matrix of dimensions
!                     (2 * (nmax + 1) + 1, 2 * (nmax + 1)) with the following
!                     structure
!
!                     [t(lativ(1),loniv(1))            t(lativ(1),loniv(2))            ...            t(lativ(1),loniv(2*(nmax+1)))]
!                     [t(lativ(2),loniv(1))            t(lativ(2),loniv(2))            ...            t(lativ(2),loniv(2*(nmax+1)))]
!                     [       .                                                    .                         .                     ]
!                     [       .                                                     .                        .                     ]
!                     [       .                                                      .                       .                     ]
!                     [t(lativ(2*(nmax+1)+1),loniv(1)) t(lativ(2*(nmax+1)+1),loniv(2)) ... t(lativ(2*(nmax+1)+1),loniv(2*(nmax+1)))]
!
!                     where "lativ" and "loniv" denote the output vectors from
!                     the "DH_grid" function.
!
!         "radius" -- Radius of the reference sphere (metres) at which the data
!                     are sampled and SRBFs are formally located; the output
!                     expansion coefficients "a" refer to this sphere; a scalar
!
!         "nmax"   -- Maximum harmonic degree up to which the analysis is
!                     performed, nmax >= 0; a scalar
!
!
! OUTPUTS: "lati" -- Spherical latitudes in radians of the nodal points, to
!                    which expansion coefficients refer to; dimensions of
!                    the vector ((2 * (nmax + 1) + 1) * (2 * (nmax + 1)), 1)
!
!          "loni" -- Spherical longitudes in radians of the nodal points, to
!                    which expansion coefficients refer to; dimensions of
!                    the vector (( 2 * (nmax + 1) + 1) * (2 * (nmax + 1)), 1)
!
!          "a"    -- Expansion coefficients (units are equal to the units
!                    of the input data but multiplied by m^2 in addition)
!                    at the nodal points (lati, loni, ri); dimensions of the
!                    vector ((2 * (nmax + 1) + 1) * (2 * (nmax + 1)), 1)
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

    integer,  intent(in)  :: nmax
    real(dp), intent(in)  :: lativ(2 * (nmax + 1) + 1)
    real(dp), intent(in)  :: loniv(2 * (nmax + 1))
    real(dp), intent(in)  :: t(2 * (nmax + 1) + 1, 2 * (nmax + 1))
    real(dp), intent(in)  :: radius

    real(dp), intent(out) :: a((2 * (nmax + 1) + 1) * (2 * (nmax + 1)))
    real(dp), intent(out) :: lati((2 * (nmax + 1) + 1) * (2 * (nmax + 1)))
    real(dp), intent(out) :: loni((2 * (nmax + 1) + 1) * (2 * (nmax + 1)))

    integer  :: i, j, l, k, n, m
    real(dp) :: c, d, colatii, w
    ! ===============================================


    ! The main computational part of the routine
    ! ===============================================
    n = 2 * (nmax + 1) + 1
    m = n - 1

    l = nmax + 1

    c = 2.0_dp * pi * radius * radius / (dble(l) * dble(l))

    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(colatii, w, d, i, k)
    ! The following "do loop" is parallelized using the number of cores
    ! specified in the "Test_run.f95" using the "maxthreads" variable. The best
    ! improvement due to the parallelization can be achieved for high values of
    ! "nmax", while for small "nmax" values, the performance may be even
    ! slightly worse due to the overhead issue. However, in the latter case,
    ! the overall walltime is only a few seconds or even less, so the increased
    ! time can be considered as negligible. Hence, the loop is parallelized
    ! by default.
    do i = 1, n

        colatii = pi / 2.0_dp - lativ(i)

        w = 0.0_dp
        do k = 0, (l - 1)
            d = (2.0_dp * dble(k) + 1.0_dp)

            w = w + c * sin(colatii) / d * sin(d * colatii)
        end do

        do j = 1, m
            a(i + (j - 1) * n)    = w * t(i, j)
            lati(i + (j - 1) * n) = lativ(i)
            loni(i + (j - 1) * n) = loniv(j)
        end do

    end do
    !$OMP END PARALLEL DO
    ! ===============================================

end subroutine SRBFs_analysis
