program Test_run
!
! ===================================================
!
! This is a Fortran program to test the SRBF_package. The program consists of
! two examples that demonstrate the correctness of the package
! in a closed-loop environment.
!
! The loops inside the "SRBFs_analysis.f95" and "SRBFs_synthesis.f95" routines
! are paralellized using OpenMP.
!
! NOTE: The user has to manually set the path to the directory, where the
!       attached input data files (*.txt) are stored (see the variable "path"
!       below).
!
! ===================================================


    ! Import modules
    ! ===============================================
    use vartypes   ! Module defining floating point double precision numbers
    use constants  ! Module defining constants (pi, etc.)
    !$ use omp_lib, only: omp_get_max_threads ! Module to get the maximum
                                              ! number of available threads. By
                                              ! default, the maximum number is
                                              ! used for parallel parts of the
                                              ! package. However, the value can
                                              ! be changed in the code below 
                                              ! (look for 
                                              ! "omp_set_num_threads").
    ! ===============================================


    ! Declaration of variables
    ! ===============================================
    implicit none



    integer :: nmin,  nmax  ! Minimum and maximum harmonic degrees of
                            ! the expansion in the EXAMPLE NO. 1, respectively
    integer :: nmin1, nmax1 ! Minimum and maximum harmonic degrees of
                            ! the expansion in the EXAMPLE NO. 2, respectively

    real(dp), allocatable :: lati(:) ! Latitudes of the nodal points, at which
                                     ! spherical radial basis functions are
                                     ! located
    real(dp), allocatable :: loni(:) ! Longitudes of the nodal points, at which
                                     ! spherical radial basis functions are
                                     ! located
    real(dp)              :: ri      ! Spherical radius of nodal points,
                                     ! constant for all spherical radial basis
                                     ! functions

    integer               :: npi     ! Number of nodal points
    integer               :: nlati   ! Number of latitudes of computational
                                     ! points at a single meridian
    integer               :: nloni   ! Number of longitudes of computational
                                     ! points at a single parallel

    real(dp), allocatable :: lat(:), lon(:), r(:) ! Latitudes, longitudes and
                                                  ! radii of all computational
                                                  ! points, respectively
    integer               :: np                   ! Number of evaluation points

    real(dp), allocatable :: a(:)       ! Expansion coefficients of the input
                                        ! signal
    integer               :: nphisrbf   ! Length of the vector "phisrbf"
    real(dp), allocatable :: phisrbf(:) ! Shape coefficients of spherical
                                        ! radial basis functions
    real(dp), allocatable :: tin(:,:)   ! Input signal from which expansion
                                        ! coefficients "a" are estimated in
                                        ! EXAMPLE NO. 1

    real(dp) :: radius ! Radius of the reference sphere, at which the nodal
                       ! points are placed
    real(dp) :: psimax ! Integration radius

    real(dp), allocatable :: lativ(:) ! Longitudes of computational points
                                        ! at a single meridian
    real(dp), allocatable :: loniv(:) ! Longitudes of computational points
                                        ! at a single meridian

    integer :: fid, i, j  ! Other useful variables

    real(dp) :: mean

    real(dp), allocatable :: txyref(:), tzzref(:) ! Reference values for the
                                                  ! EXAMPLE NO. 1 and 2,
                                                  ! respectively

    real(dp), allocatable :: diffs(:)  ! Differences to compute the statistics

    character(len = 256) :: path ! Absolute path to the directory, where data
                                 ! files
                                 !
                                 ! (i)   "Input_Data.txt",
                                 ! (ii)  "Example1_ref_data.txt", and
                                 ! (iii) "Example2_ref_data.txt"
                                 !
                                 ! are stored.
                                 !
                                 ! This variable needs to be set correctly
                                 ! by the user to run this test program (see
                                 ! below).

    integer :: nmaxthreads ! Maximum number of threads used for parallel parts
                           ! of the package ("SRBFs_analysis.f95" and
                           ! "SRBFs_synthesis.f95")

    real(dp), allocatable ::   t(:)                 ! Disturbing potential
    real(dp), allocatable ::  tx(:),  ty(:),  tz(:) ! Disturbing vector
    real(dp), allocatable :: txx(:), txy(:), txz(:) ! Disturbing tensor
    real(dp), allocatable :: tyy(:), tyz(:), tzz(:) ! Disturbing tensor
    ! ===============================================


    ! The main computational part of the program
    ! ===============================================
    ! Absolute path to the directory, where data files
    !
    ! (i)   "Input_Data.txt",
    ! (ii)  "Example1_ref_data.txt", and
    ! (iii) "Example2_ref_data.txt"
    !
    ! are stored.
    !
    ! The "path" variable needs to be set correctly by the user to run this test
    ! program.
    path = '<Root directory to SRBF_package>/SRBF_package/Fortran/Input_data'

    nmin = 0
    nmax = 35

    radius = 6378136.3_dp
    ri     = radius
    psimax = pi

    nlati = (2 * (nmax + 1) + 1)
    nloni = nlati - 1
    npi   = nlati * nloni

    ! Get and then set the maximum number of available threads for parallel
    ! parts of the package.
    !$ nmaxthreads = omp_get_max_threads()
    !$ call omp_set_num_threads(nmaxthreads)




    ! TEST EXAMPLE NO. 1
    ! In this example, the input quantity is the disturbing potential
    ! synthesized from global geopotential model EGM96 via spherical harmonic
    ! synthesis (SHS) in MATLAB. The minimum and the maximum degrees of the
    ! synthesis were "nmin" and "nmax" (see above), respectively. The disturbing
    ! potential is synthesized on the reference sphere with the radius "radius".
    ! These are the input data for SRBFs analysis. Through the SRBFs analysis,
    ! the expansion coefficients "a" are estimated via the exact quadrature.
    ! Finally, the Txy element of the disturbing tensor is synthesized
    ! from the expansion coefficients (via SRBF synthesis). To check the
    ! correctness of the upward continuation operations, the validation is
    ! performed on a sphere being 1000 m above the reference sphere. These
    ! results are then compared with reference values from SHS. The differences
    ! between the two results reveal the associated errors of the SRBFs
    ! analysis and synthesis (errors in the reference data from SHS are here
    ! neglected).
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    print *, '==================================='
    print *, 'TEST EXAMPLE NO. 1'
    print *, ''

    ! Read input data from the "Input_Data.txt" file.
    ! -----------------------------------------------
    ! This file contains the disturbing potential from EGM96 obtained by the
    ! "GrafLab.m" software using the following parameters (the same ones as
    ! used in the MATLAB version of the SRBF_package):
    !
    ! nmin = 0
    ! nmax = 35
    ! R = 6378136.3 (radius of the reference sphere related to EGM96)
    ! GM = 3986004.415E+8 (geocentric gravitational constant related to EGM96)
    ! lativ, loniv (spherical latitudes and longitudes of the evaluation
    !                   points that can be obtained by as
    !                   "DH_grid(35, lativ, loniv)", see below)
    ! r = R (radius of the evaluation points)
    ! GRS80 elipsoid to generate the normal gravitational field

    print *, 'Loading input data...'

    open(newunit = fid, file = trim(adjustl(path))//'/Input_data.txt', &
         status = 'old')

    allocate(tin(nlati, nloni))

    do i = 1, nlati
        read(fid, *) tin(i, :)
    end do
    ! -----------------------------------------------


    ! Compute the latitudes "lativ" and longitudes "loniv" of the grid,
    ! at which input data from "Input_data.txt" file are given.
    ! -----------------------------------------------
    print *, 'Computing the Driscoll--Healy grid...'

    allocate(lativ(nlati), loniv(nloni))

    call DH_grid(nmax, &
                 lativ, loniv)
    ! -----------------------------------------------


    ! Now compute the latitudes "lat" and longitudes "lon" of the evaluation
    ! points, at which the validation of the package is to be performed. Here,
    ! we simply use the same grid grid boundaries as defined by the nodal points
    ! "lativ" and "loniv". A regular grid is created from these boudnaries.
    ! All evaluation points share the same constant radius of "radius + 1000 m".
    ! -----------------------------------------------
    print *, 'Creating grid of evaluation points...'

    np = npi ! Number of evaluation points is equal to the total number of
             ! nodal points (in general, can take any value that is equal
             ! to or larger than 1)

    allocate(lat(np), lon(np), r(np))

    do i = 1, nlati
        do j = 1, nloni
            lat(i + (j - 1) * nlati) = lativ(i)
            lon(i + (j - 1) * nlati) = loniv(j)
        end do
    end do

    r = radius + 1000.0_dp
    ! -----------------------------------------------


    ! Perform the SRBFs analysis
    ! -----------------------------------------------
    print *, 'SRBFs analysis...'

    allocate(a(npi), lati(npi), loni(npi))

    call SRBFs_analysis(lativ, loniv, tin, radius, nmax, &
                        lati, loni, a)
    ! -----------------------------------------------


    ! Perform the SRBFs synthesis
    ! -----------------------------------------------
    print *, 'SRBFs synthesis...'

    nphisrbf = nmax + 1

    allocate(phisrbf(nphisrbf), &
             t(np), tx(np), ty(np), tz(np),&
             txx(np), txy(np), txz(np),  &
             tyy(np), tyz(np), tzz(np))

    phisrbf = 1.0_dp ! Shannon spherical radial basis function is used in
                     ! the EXAMPLE NO. 1

    call SRBFs_synthesis(lati, loni, ri, npi, lat, lon, r, np, nmin, nmax, a, &
                         phisrbf, nphisrbf, radius, psimax, &
                         t, tx, ty, tz, txx, txy, txz, tyy, tyz, tzz)
    ! -----------------------------------------------


    ! Read reference data for the element txy that serve as sample reference
    ! data for the validation
    ! -----------------------------------------------
    print *, 'Loading the reference data (the Txy element)'

    open(newunit = fid, file = trim(adjustl(path))//'/Example1_ref_data.txt', &
         status = 'old')

    allocate(txyref(np))

    do i = 1, np
        read(fid, *) txyref(i)
    end do
    ! -----------------------------------------------


    ! Perform the validation
    ! -----------------------------------------------
    allocate(diffs(np))

    ! Difference between the SRBF synthesis and reference values
    diffs = (txy * 1e9_dp) - txyref

    print *, ''
    print *, 'Statistics of the differences between the SRBF synthesis and &
              &the references values'

    print *, 'Min. (E)'
    write(*,*) minval(diffs)

    print *, ''
    print *, 'Max. (E)'
    write(*,*) maxval(diffs)

    print *, ''
    print *, 'Mean (E)'
    mean = sum(diffs) / dble(size(diffs))
    write(*,*) mean

    print *, ''
    print *, 'STD (E)'
    write(*,*) sqrt(sum((diffs - mean)**2) / (size(diffs)))

    print *, ''
    print *, '(These values should be generally below the 10**(-13) E level)'
    ! -----------------------------------------------
    print *, '==================================='
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! End of TEST EXAMPLE NO. 1






    ! TEST EXAMPLE NO. 2
    ! In this example, the expansion coefficients from "Test Example No. 1" are
    ! used to synthesize the Tzz element of the disturbing tensor, but this time
    ! within a different spectral band than that for which the coefficients
    ! where estimated. More specifically, while the coefficients cover the
    ! harmonic degrees "nmin" -- "nmax", here, we use them to synthesize Tzz
    ! within the spectral band "nmin1" -- "nmax1", where "nmin < nmin1 <= nmax"
    ! and "nmin1 <= nmax1 <= nmax". The Tzz values synthesized from SRBFs are
    ! then compared with SHS as a reference. The validation is conducted on a
    ! sphere being 1000 m above the reference sphere.
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    print *, ''
    print *, ''
    print *, ''
    print *, '==================================='
    print *, 'TEST EXAMPLE NO. 2'
    print *, ''

    nmin1 = 10 ! Minumum degree of the synthesis
    nmax1 = 20 ! Maximum degree of the synthesis

    deallocate(phisrbf)
    nphisrbf = nmax1 + 1
    allocate(phisrbf(nphisrbf))

    phisrbf = 1.0_dp ! Shannon spherical radial basis function is used
                     ! in the EXAMPLE NO. 2

    ! Perform the SRBFs synthesis
    ! -----------------------------------------------
    print *, 'SRBFs synthesis...'

    call SRBFs_synthesis(lati, loni, ri, npi, lat, lon, r, np, nmin1, nmax1, &
                         a, phisrbf, nphisrbf, radius, psimax, &
                         t, tx, ty, tz, txx, txy, txz, tyy, tyz, tzz)
    ! -----------------------------------------------


    ! Read reference data for the element tzz that serve as sample reference
    ! data for the validation
    ! -----------------------------------------------
    print *, 'Loading the reference data (the Tzz element)'

    open(newunit = fid, file = trim(adjustl(path))//'/Example2_ref_data.txt', &
         status = 'old')

    allocate(tzzref(np))

    do i = 1, np
        read(fid, *) tzzref(i)
    end do
    ! -----------------------------------------------


    ! Perform the validation
    ! -----------------------------------------------

    ! Difference between the SRBF synthesis and reference values
    diffs = (tzz * 1e9_dp) - tzzref

    print *, ''
    print *, 'Statistics of the differences between the SRBF synthesis and &
              &the references values'
    print *, 'Min. (E)'
    write(*,*) minval(diffs)

    print *, ''
    print *, 'Max. (E)'
    write(*,*) maxval(diffs)

    print *, ''
    print *, 'Mean (E)'
    mean = sum(diffs) / dble(size(diffs))
    write(*,*) mean

    print *, ''
    print *, 'STD (E)'
    write(*,*) sqrt(sum((diffs - mean)**2) / (size(diffs)))
    ! -----------------------------------------------

    print *, ''
    print *, '(These values should be generally below the 10**(-12) E level)'
    print *, '==================================='
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! End of TEST EXAMPLE NO. 2




    print *, ''
    print *, ''
    print *, ''
    print *, 'End of the test run.'

end program Test_run
