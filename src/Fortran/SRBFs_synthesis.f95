subroutine SRBFs_synthesis(lati, loni, ri, npi, lat, lon, r, np, nmin, nmax, &
                           a, phisrbf, nphisrbf, radius, psimax, &
                           v, vx, vy, vz, vxx, vxy, vxz, vyy, vyz, vzz)
!
! =========================================================================
!
! DESCRIPTION: This function synthesizes the gravitational potential, the
!              gravitational vector and the gravitational tensor from a given
!              set of expansion coefficients related to band-limited spherical
!              radial basis functions (SRBFs).
!
!              For the adopted definition of SRBFs, see the attached PDF
!              "Gravitational_potential_up_to_its_second-order_derivatives_in_
!              terms_of_SRBFs.pdf".
!              The gravitational vector and tensor are expressed in the local
!              north-oriented reference frame (LNOF). The definition of LNOF
!              as well as of the formulae used for the synthesis are provided
!              in the attached PDF.
!
!              The loop over evaluation points is parallelized using OpenMP.
!
! INPUTS: "lati"     -- Spherical latitudes in radians of nodal points, at which
!                       SRBFs are formally located; a vector of dimensions
!                        (npi, 1)
!
!         "loni"     -- Spherical longitudes in radians of nodal points, at
!                       which SRBFs are formally located; a vector of dimensions
!                       (npi, 1)
!
!         "ri"       -- Spherical radius (in metres) of the sphere, at which
!                       SRBFs are formally located; a scalar ("ri" is assumed
!                       to be constant for all nodal points)
!
!         "npi"      -- Number of nodal points
!
!         "lat"      -- Spherical latitudes (in radians) of evaluation points;
!                       a vector of dimensions (np, 1)
!
!         "lon"      -- Spherical longitudes (in radians) of evaluation points;
!                       a vector of dimensions (np, 1)
!
!         "r"        -- Spherical radii in metres of evaluation points;
!                       a vector of dimensions (np, 1)
!
!         "np"       -- Number of evaluation points
!
!         "nmin"     -- Minimum harmonic degree of the synthesis; a scalar
!                       satisfying the inequality 0 <= nmin <= nmax
!
!         "nmax"     -- Maximum harmonic degree of the synthesis; a non-negative
!                       scalar
!
!         "a"        -- Expansion coefficients of SRBFs in m^4 * s^-2 (in the
!                       case of the input gravitational potential
!                       in m^2 * s^-2); a vector of dimensions (npi, 1)
!
!         "phisrbf"  -- Shape coefficients of SRBFs (dimensionless); a vector of
!                       dimensions (nphisrbf, 1)
!
!         "nphisrbf" -- Length of the vector "phisrbf"
!
!         "radius"   -- Radius of the reference sphere, at which the nodal
!                       points (ri, lati, loni) are placed (metres); a scalar
!
!         "psimax"   -- Spherical distance in radians, up to which SRBFs
!                       placed at the nodal points (lati, loni, ri) are
!                       considered in the synthesis; beyond psimax, the
!                       expansion coefficients of SRBFs at the nodal points
!                       (lati, loni, ri) are automatically assumed to be equal
!                       to zero; for global integration, choose "psimax = pi",
!                       for regional (cap) integration, choose
!                       "0 < psimax < pi"; a scalar
!
!
! OUTPUTS: "v"   -- Gravitational potential synthesized at the evaluation
!                   points (lat, lon, r) in m^2 * s^-2; a vector of dimensions
!                   (np, 1)
!
!          "vx, vy, vz" -- Elements of the gravitational vector in LNOF
!                   synthesized at the evaluation points (lat, lon, r) in
!                   m * s^-2; three vectors of dimensions (np, 1)
!
!          "vxx, vxy, vxz, vyy, vyz, vzz" -- Elements of the gravitational
!                   tensor in LNOF synthesized at the evaluation points
!                   (lat, lon, r) in s^-2; six vectors of dimensions (np, 1)
!
!
! REFERENCES: Bucha, B., Bezdek, A., Sebera, J., Janak, J., 2015. Global and
!                regional gravity field determination from GOCE kinematic orbit
!                by means of spherical radial basis functions. Surveys in
!                Geophysics 36, 773-801,
!                http://doi.org/10.1007/s10712-015-9344-0
!                (preprint freely available at
!                https://www.researchgate.net/profile/Blazej_Bucha)
!
!             Bucha, B., Janak, J., Papco, J., Bezdek, A., 2016. High-resolution
!                regional gravity field modelling in a mountainous area
!                from terrestrial gravity data. Geophysical Journal
!                International 207, 949-966, http://doi.org/10.1093/gji/ggw311
!                (preprint freely available at
!                https://www.researchgate.net/profile/Blazej_Bucha)
!
!             Bucha, B., Hirt, C., Kuhn, M., 2019. Cap integration in spectral
!                gravity forward modelling up to the full gravity tensor.
!                Journal of Geodesy, https://doi.org/10.1007/s00190-019-01277-3
!                (preprint freely available at
!                https://www.researchgate.net/profile/Blazej_Bucha)
!
!             Sprlak, M., Hamackova, E., Novak, P., 2015 Alternative
!                validation method of satellite gradiometric data by integral
!                transform of satellite altimetry data. Journal of Geodesy 89,
!                757-773, http://doi.org/10.1007/s00190-015-0813-5
!
!
! Please use the following reference when using this function:
!
!    Bucha, B., Janak, J., Papco, J., Bezdek, A., 2016. High-resolution
!        regional gravity field modelling in a mountainous area
!        from terrestrial gravity data. Geophysical Journal
!        International 207, 949-966, http://doi.org/10.1093/gji/ggw311
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


    ! Parameters of the expansion
    ! ---------------------------------
    integer,  intent(in) :: nmin    ! Minimum degree of the expansion
    integer,  intent(in) :: nmax    ! Maximum degree of the expansion
    real(dp), intent(in) :: radius  ! Radius of the reference sphere
    ! ---------------------------------


    ! Nodal points
    ! ---------------------------------
    integer,  intent(in) :: npi                       ! Number of nodal points
    real(dp), intent(in) :: lati(npi), loni(npi)      ! latitude (rad),
                                                      ! longitude (rad),
                                                      ! respectively
    real(dp), intent(in) :: ri                        ! radius (m), constant for
                                                      ! all nodal points
    real(dp)             :: xi(npi), yi(npi), zi(npi) ! cartesian coordinates of
                                                      ! nodal points residing on
                                                      ! the unit sphere
    real(dp), intent(in) :: psimax                    ! Integration radius in
                                                      ! radians (distance up to
                                                      ! which nodal points are
                                                      ! used in the synthesis)
    ! ---------------------------------


    ! Evaluation points
    ! ---------------------------------
    integer,  intent(in) :: np                      ! Number of evaluation
                                                    ! points, respectively
    real(dp), intent(in) :: lat(np), lon(np), r(np) ! Evaluation points:
                                                    ! latitude (rad),
                                                    ! longitude (rad),
                                                    ! radius (m), respectively
    real(dp)             :: x(np), y(np), z(np)     ! Evaluation points:
                                                    ! cartesian coordinates
    ! ---------------------------------


    ! Spherical radial basis functions
    ! ---------------------------------
    real(dp), intent(in) :: a(npi)            ! Expansion coefficients
    real(dp)             :: ai                ! The ith expansion coefficient

    integer,  intent(in) :: nphisrbf             ! Length of the vector
                                                 ! "phisrbf"
    real(dp), intent(in) :: phisrbf(nphisrbf)    ! Shape coefficients of
                                                 ! spherical radial basis
                                                 ! functions
    real(dp)             :: phisrbftmp(nphisrbf) ! Temporary "phisrbf"

    real(dp) :: phi,   phi10, phi11 ! Spherical radial basis functions Phi,
                                    ! Phi10, Phi11
    real(dp) :: phi20, phi21, phi22 ! Spherical radial basis functions Phi20,
                                    ! Phi21, Phi22

    ! Substitutions for coefficients of degree n of spherical radial
    ! basis functions
    real(dp) ::   phin(nmax + 1), phin10(nmax + 1), phin11(nmax + 1)
    real(dp) :: phin20(nmax + 1), phin21(nmax + 1), phin22(nmax + 1)
    ! ---------------------------------


    ! Useful substitutions
    ! ---------------------------------
    real(dp) :: cpsi, spsi       ! cos and sin of the spherical distance between
                                 ! the jth computation point and all nodal
                                 ! points, respectively
    real(dp) :: spsi2            ! spsi**2
    real(dp) :: alpha            ! Azimuth between the jth evaluation point and
                                 ! all nodal points
    real(dp) :: calpha, salpha   ! cos(alpha), sin(alpha)
    real(dp) :: c2alpha, s2alpha ! cos(2 * alpha), sin(2 * alpha)
    real(dp) :: t, tpow2, tpown  ! t = ri / r(j), t**2 and t**n, respetively


    real(dp) :: c(nmax + 1), d(nmax + 1)  ! Coefficients to compute Legendre
    real(dp) :: e(nmax + 1), f(nmax +  1) ! polynomials
    real(dp) :: p0, p1, p2     ! Legendre polynomials of degrees "n - 2",
                               ! "n - 1" and "n", respectively
    real(dp) :: dp1, dp2       ! First-order derivatives of Legendre polynomials
    real(dp) :: ddp1, ddp2     ! Second-order derivatives of Legendre
                               ! polynomials
    real(dp) :: vxx1, vxx2     ! Auxiliary variable to compute the Vxx and Vyy
                               ! elements of the gravitational tensor

    real(dp) :: r2pi4   ! 4.0_dp * pi * radius**2
    real(dp) :: radius2 ! radius**2

    integer  :: i, j, n ! Loop variables for nodal points, evaluation points
                        ! and harmonic degrees, respectively

    real(dp) :: cpsimax
    ! ---------------------------------


    ! Output variables
    ! ---------------------------------
    ! Gravitational potential
    real(dp), intent(out) ::   v(np)

    ! Gravitational vector
    real(dp), intent(out) ::  vx(np),  vy(np),  vz(np)

    ! Gravitational tensor
    real(dp), intent(out) :: vxx(np), vxy(np), vxz(np)
    real(dp), intent(out) :: vyy(np), vyz(np), vzz(np)
    ! ---------------------------------
    ! ===============================================


    ! The main computational part of the routine
    ! ===============================================

    ! Initialization of output variables
    ! ----------------------------------------------
    v   = 0.0_dp
    vx  = 0.0_dp
    vy  = 0.0_dp
    vz  = 0.0_dp
    vxx = 0.0_dp
    vxy = 0.0_dp
    vxz = 0.0_dp
    vyy = 0.0_dp
    vyz = 0.0_dp
    vzz = 0.0_dp
    ! ----------------------------------------------


    ! Transformation of spherical coordinates of the nodal and evaluation
    ! points to cartesian coordinates, assuming the unit sphere
    ! ----------------------------------------------
    call sph2cart(loni, lati, npi, xi, yi, zi)
    call sph2cart(lon,  lat,  np,   x,  y,  z)
    ! ----------------------------------------------


    ! Computation of coefficients for recurrence relations related to Legendre
    ! polynomials (LP) and their derivatives
    ! ----------------------------------------------
    c(1) = 0.0_dp ! Non-existing coefficient set to zero
    d(1) = 0.0_dp ! Non-existing coefficient set to zero
    e(1) = 0.0_dp ! Existing coefficient, but is equal to zero
    f(1) = 1.0_dp
    do n = 1, nmax

        c(n + 1) = (2.0_dp * dble(n) - 1.0_dp) / dble(n)
        d(n + 1) = (dble(n) - 1.0_dp) / dble(n)
        e(n + 1) = dble(n)
        f(n + 1) = dble(n) + 1.0_dp

    end do
    ! ----------------------------------------------


    ! Substitutions
    ! ----------------------------------------------
    radius2 = radius * radius
    r2pi4   = 4.0_dp * pi * radius2
    cpsimax = cos(psimax)
    ! ----------------------------------------------


    phisrbftmp = phisrbf
    ! If "nmin > 0", the shape coefficients of SRBFs are set to zero
    ! to ensure correct SRBF synthesis
    phisrbftmp(1:nmin) = 0.0_dp


    ! Loop over evaluation points
    ! ----------------------------------------------
    ! The following "do loop" is parallelized using the number of cores
    ! specified in the "Test_run.f95" using the "maxthreads" variable.
    !$OMP PARALLEL DO DEFAULT(PRIVATE) &
    !$OMP& SHARED(lat, lon, r, x, y, z, np, lati, loni, ri, xi, yi, zi, npi) &
    !$OMP& SHARED(nmax, psimax, cpsimax, radius, radius2, r2pi4, phisrbftmp) &
    !$OMP& SHARED(a, c, d, e, f, v, vx, vy, vz, vxx, vxy, vxz, vyy, vyz, vzz)
    do j = 1, np

        t = ri / r(j)
        tpow2 = t * t
        tpown = t

        ! Loop over harmonic degrees to prepare coefficients of spherical
        ! radial basis functions
        do n = 0, nmax
            ! Substitutions for radial basis function Phi
            phin(n + 1) = (2.0_dp * dble(n) + 1.0_dp) / r2pi4 * &
                          tpown * phisrbftmp(n + 1)
            tpown = tpown * t


            ! Substitutions for radial basis function Phi10 and Phi11
            phin10(n + 1) = -1.0_dp / radius * (dble(n) + 1.0_dp) * &
                             t * phin(n + 1)
            phin11(n + 1) = -1.0_dp / radius * t * phin(n + 1)


            ! Substitutions for radial basis function Phi20, Phi21 and Phi22

            phin20(n + 1) =  1.0_dp / radius2 * (dble(n) + 1.0_dp) * &
                             (dble(n) + 2.0_dp) * tpow2 * phin(n + 1)
            phin21(n + 1) = -1.0_dp / radius2 * (dble(n) + 2.0_dp) * &
                             tpow2 * phin(n + 1)
            phin22(n + 1) =  1.0_dp / (2.0_dp * radius2) * tpow2 * phin(n + 1)
        end do


        do i = 1, npi

            ! Spherical distance psi and sin(psi), cos(psi), etc.
            ! -------------------------------------------------
            cpsi = x(j) * xi(i) + y(j) * yi(i) + z(j) * zi(i) ! cos(psi)

            if (abs(pi - psimax) >= 1e-14_dp) then ! Regional integration
                if (cpsi < cpsimax) then ! If the nodal point is beyond
                                         ! the distance specified via
                                         ! the "psimax" variable,
                                         ! the contribution from that nodal
                                         ! point is omitted
                    cycle
                end if
            end if

            ! To avoid numerical inaccuracies that might cause
            ! imaginary numbers in the output quantities
            if (cpsi > 1.0_dp) then
                cpsi = 1.0_dp
            elseif (cpsi < -1.0_dp) then
                cpsi = -1.0_dp
            end if

            spsi  = sin(acos(cpsi)) ! sin(psi)
            spsi2 = spsi * spsi     ! sin(psi2)**2
            ! ---------------------------------------------


            ! Azimuth "alpha"
            ! ---------------------------------------------
            alpha   = atan2(cos(lati(i)) * sin(loni(i)  - lon(j)), &
                            cos(lat(j))  * sin(lati(i)) - sin(lat(j)) * &
                            cos(lati(i)) * cos(loni(i)  - lon(j)))
            calpha  = cos(alpha)
            salpha  = sin(alpha)
            c2alpha = calpha * calpha - salpha * salpha ! sin(2.0_dp * alpha(i))
            s2alpha = 2.0_dp * salpha * calpha          ! cos(2.0_dp * alpha(i))
            ! ---------------------------------------------


            ! Computation of SRBFs Phi, Phi10, Phi11, Phi20, Phi21, Phi22
            ! ---------------------------------------------
            ! Seed values for Legendre polynomials
            ! .............................................
            p0   = 1.0_dp
            p1   = cpsi
            dp1  = p0
            ddp1 = 0.0_dp
            ! .............................................


            ! Sum over harmonic degrees 0 and 1
            ! .............................................
            phi   = phin(1) * p0
            phi10 = phin10(1) * p0
            phi20 = phin20(1) * p0

            if (nmax > 0) then
                phi   = phi + phin(2) * p1
                phi10 = phi10 + phin10(2) * p1
                phi11 = phin11(2) * dp1
                phi20 = phi20 + phin20(2) * p1
                phi21 = phin21(2) * dp1
                phi22 = 0.0_dp
            else
                phi11 = 0.0_dp
                phi21 = 0.0_dp
                phi22 = 0.0_dp
            end if
            ! .............................................


            ! Sum over harmonic degrees 2..nmax
            ! .............................................
            do n = 2, nmax

                ! Recurrence relation for LPs
                p2  = c(n + 1) * cpsi * p1 - d(n + 1) * p0

                ! Recurrence relation for dLPs
                dp2 = e(n + 1) * p1 + cpsi * dp1

                ! Recurrence relation for ddLPs
                ddp2 = f(n + 1) * dp1 + cpsi * ddp1

                phi   = phi   +   phin(n + 1) *   p2
                phi10 = phi10 + phin10(n + 1) *   p2
                phi11 = phi11 + phin11(n + 1) *  dp2
                phi20 = phi20 + phin20(n + 1) *   p2
                phi21 = phi21 + phin21(n + 1) *  dp2
                phi22 = phi22 + phin22(n + 1) * ddp2

                p0   =   p1
                p1   =   p2
                dp1  =  dp2
                ddp1 = ddp2

            end do ! End of the loop over harmonic degrees
            ! .............................................


            ! Computation of radial basis functions Phi11, Phi21 and Phi22
            ! .............................................
            phi11 = phi11 * spsi
            phi21 = phi21 * spsi
            phi22 = phi22 * spsi2
            ! .............................................


            ! Synthesis with spherical radial basis functions
            ! .............................................
            ai = a(i)

            vxx1 = -0.5_dp * ai * phi20
            vxx2 = ai * c2alpha * phi22

            v(j)   = v(j)   +   ai * phi
            vx(j)  = vx(j)  + (-ai * phi11 * calpha)
            vy(j)  = vy(j)  +   ai * phi11 * salpha
            vz(j)  = vz(j)  +   ai * phi10
            vxx(j) = vxx(j) + vxx1  + vxx2
            vxy(j) = vxy(j) + (-ai * phi22 * s2alpha)
            vxz(j) = vxz(j) +   ai * phi21 * calpha
            vyy(j) = vyy(j) + vxx1 - vxx2
            vyz(j) = vyz(j) + (-ai * phi21 * salpha)
            vzz(j) = vzz(j) +  (ai * phi20)
            ! .............................................

        end do

    end do ! End of loop over the evaluation points
    !$OMP END PARALLEL DO
    ! ----------------------------------------------
    ! ===============================================

end subroutine SRBFs_synthesis
