# -*- coding: utf-8 -*-
"""
===================================================

This is a Python script to test the Fortran routines of the SRBF_package
in Python. The script consists of two examples that demonstrate the correctness
of the package in a closed-loop environment.

Before running the script, three remarks are necessary.

1) In order to be able to employ the Fortran routines in Python, a Python
   wrapper has to be created. An easy way to do this is to use F2PY, which is
   a part of the Python's "numpy" module. For

   1) Debian GNU/Linux operating system,
   2) GNU Fortran compiler (Debian 8.3.0-6) 8.3.0,
   3) OpenMP,
   4) Python 3.7.3,
   5) numpy 1.19.0.dev0+f71d993,

   it was tested that the wrapper can correctly be created by the following
   command (assiming that your current working directory is
   <Root_Directory_of_the_Package/SRBF_package/Fortran/Python>):

   f2py -c ../vartypes.f95 ../constants.f95 ../SRBFs_synthesis.f95 \
   ../SRBFs_analysis.f95 ../DH_grid.f95 ../sph2cart.f95 -m SRBF_package \
   --f90flags='-fopenmp -O3 -ffast-math -expensive-optimizations -funroll-loops' \
   -lgomp

   As a result, a wrapper named "SRBF_package" will be created in the
   "SRBF_package/Fortran/Python/" folder. After successfully creating
   the wrapper, the "Test_run.py" script can be executed.

   The signature file "SRBF_package.pyf", from which the wrapper can
   alternatively be obtained, is available in "SRBF_package/Fortran/Python".

   Two other options to prepare the wrapper can be found at
   https://www.fortran90.org/src/best-practices.html#interfacing-with-python

2) "do loops" inside the "SRBFs_analysis.f95" and "SRBFs_synthesis.f95"
   Fortran subroutines are paralellized using OpenMP. By default, the maximum
   number of available threads is used. If some other value is required,
   it can be specified by the environment variable "OMP_NUM_THREADS"
   by the user.

3) Info on the inputs and outputs of the functions from the "SRBF_package"
   module can be found using commands such as

   print(SP.srbfs_synthesis.__doc__)

   which prints the inputs and outputs of the "srbfs_synthesis" function from
   the "SRBF_package" module, etc. For most of the functions from the module,
   the same syntax as with the MATLAB counterpart of the toolbox can be used.
   However, in Python some inputs may be additional and/or optional.

   Further details on the inputs and outputs can also be found in the signature
   file "SRBF_package.pyf", which is available in the
   "SRBF_package/Fortran/Python" folder.

===================================================
"""




# Import modules
# ===============================================

import os
import numpy as np
import SRBF_package as SP
import matplotlib as plt

# ===============================================




# The main computational part of the routine
# ===============================================

# Set some input variables at first
# -----------------------------------------------
path = '../Input_data' # Absolute or relative path to the directory, where data
                       # files
                       #
                       # (i)   "Input_Data.txt",
                       # (ii)  "Example1_ref_data.txt", and
                       # (iii) "Example2_ref_data.txt"
                       #
                       # are stored.

nmin = 0   # Minimum harmonic degree of the expansion for EXAMPLE NO. 1
nmax = 35  # Maximum harmonic degree of the expansion for EXAMPLE NO. 1

R  = 6378136.3  # Radius of the reference sphere
ri = R          # Radius of the nodal points
pi = 3.14159265358979324
psimax = pi     # Global integration with the integration radius 180 deg
# -----------------------------------------------





# TEST EXAMPLE NO. 1
# In this example, the input quantity is the disturbing potential
# synthesized from global geopotential model EGM96 via spherical harmonic
# synthesis (SHS) in MATLAB. The minimum and the maximum degrees of the
# synthesis were "nmin" and "nmax" (see above), respectively. The disturbing
# potential is synthesized on the reference sphere with the radius "radius".
# These are the input data for SRBFs analysis. Through the SRBFs analysis,
# the expansion coefficients "a" are estimated via the exact quadrature.
# Finally, the Txy element of the disturbing tensor is synthesized
# from the expansion coefficients (via SRBF synthesis). To check the
# correctnessof the upward continuation operations, the validation is
# performed on a sphere being 1000 m above the reference sphere. These
# results are then compared with reference values from SHS. The differences
# between the two results reveal the associated errors of the SRBFs
# analysis and synthesis (errors in the reference data from SHS are here
# neglected).
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('===================================')
print('TEST EXAMPLE NO. 1\n')


# Read input data from the "Input_Data.txt" file.
# -----------------------------------------------
# This file contains the disturbing potential from EGM96 obtained by the
# "GrafLab.m" software using the following parameters (the same ones as
# used in the MATLAB version of the SRBF_package):
#
# nmin = 0
# nmax = 35
# R = 6378136.3 (radius of the reference sphere related to EGM96)
# GM = 3986004.415E+8 (geocentric gravitational constant related to EGM96)
# lativ, loniv (spherical latitudes and longitudes of the evaluation
#                   points that can be obtained by as
#                   "DH_grid(35, lativ, loniv)", see below)
# r = R (radius of the evaluation points)
# GRS80 elipsoid to generate the normal gravitational field

print('Loading input data...')

pathtmp = os.path.join(path, 'Input_data.txt')
T = np.loadtxt(pathtmp)
# -----------------------------------------------


# Compute the latitudes "lativ" and longitudes "loniv" of the grid,
# at which input data from "Input_data.txt" file are given.
# -----------------------------------------------
print('Computing the Driscoll--Healy grid...')

(lativ, loniv) = SP.dh_grid(nmax)

nlati = len(lativ)  # Number of nodal points in a single meridian
nloni = len(loniv)  # Number of nodal points in a single parallel
npi = nlati * nloni   # Number of nodal points
# -----------------------------------------------


# Plot the input data
# -----------------------------------------------
print('Plotting the input data...')

plt.pyplot.figure()
plt.pyplot.imshow(T, interpolation='none', origin='lower',
                  extent=([min(loniv * 180.0 / pi),
                           max(loniv * 180.0 / pi),
                           min(lativ * 180.0 / pi),
                           max(lativ * 180.0 / pi)]))
plt.pyplot.xlabel('Longitude (deg)')
plt.pyplot.ylabel('Latitude (deg)')
plt.pyplot.title('Example No. 1: Input quantity for the SRBFs analysis\n'
                 '(the disturbing potential in m**2 * s**(-2), '
                 'spectral band %d -- %d)' % (nmin, nmax))
plt.pyplot.colorbar(orientation='horizontal')
plt.pyplot.grid()
plt.pyplot.show()
# -----------------------------------------------


# Now compute the latitudes "lat" and longitudes "lon" of the evaluation
# points, at which the validation of the package is to be performed. Here,
# we simply use the same grid grid boundaries as defined by the nodal points
# "lativ" and "loniv". A regular grid is created from these boudnaries.
# All evaluation points share the same constant radius of "radius + 1000 m".
# -----------------------------------------------
print('Creating grid of evaluation points...')

(lon, lat) = np.meshgrid(loniv, lativ, indexing='ij')

lon = np.matrix.flatten(lon)  # Vector of latitudes of evaluation points
lat = np.matrix.flatten(lat)  # Vector of longitudes of evaluation points
r = np.zeros(lon.shape) + R + 1000.0  # Spherical radii of evaluation points

np2 = len(lon)  # Number of evaluation points
# -----------------------------------------------


# Perform the SRBFs analysis
# -----------------------------------------------
print('SRBFs analysis...')

(lati, loni, a) = SP.srbfs_analysis(lativ, loniv, T, R, nmax)
# -----------------------------------------------


# Plot the expansion coefficients
# -----------------------------------------------
print('Plotting the expansion coefficients...')

plt.pyplot.figure()
plt.pyplot.scatter(loni * 180.0 / pi, lati * 180.0 / pi, s=10, c=a)
plt.pyplot.xlabel('Longitude (deg)')
plt.pyplot.ylabel('Latitude (deg)')
plt.pyplot.title('Example No. 1: Expansion coefficients of SRBFs '
                 '(m**4 * s**(-2)')
plt.pyplot.colorbar(orientation='horizontal')
plt.pyplot.grid()
plt.pyplot.show()
# -----------------------------------------------


# Perform the SRBFs synthesis
# -----------------------------------------------
print('SRBFs synthesis...')

phisrbf = np.ones((nmax + 1, 1))  # Shape coefficients of SRBFs

(T, Tx, Ty, Tz, Txx, Txy, Txz, Tyy, Tyz, Tzz) = SP.srbfs_synthesis(lati,
 loni, ri,  lat, lon, r, nmin, nmax, a, phisrbf, R, psimax)
# Note that the syntax of "SP.srbfs_synthesis" is the same as in the MATLAB
# counterpart of the toolbox, but slightly different than in the Fortran
# subroutine "SRBFs_synthesis.f95". Here, three optional input parameters
# can be specified, defining the number of nodal points, the number of
# evaluation points and the size of the "phisrbf" vector, respectively.
#
# Using these optional parameter, the previous command reads
#
# (T, Tx, Ty, Tz, Txx, Txy, Txz, Tyy, Tyz, Tzz) = SP.srbfs_synthesis(lati,
# loni, ri,  lat, lon, r, nmin, nmax, a, phisrbf, R, psimax, npi, np2,
# len(phisrbf))
#
# and yields the same results. Further details on the inputs and outputs
# of "SP.srbfs_synthesis" can be found using the command
# print(SP.srbfs_synthesis.__doc__)
# -----------------------------------------------


# Plot the results of the synthesis
# -----------------------------------------------
print('Plotting the Txy element from the SRBF synthesis...')

Txyplt = Txy.reshape(nlati, nloni, order='F') * 1e9
plt.pyplot.figure()
plt.pyplot.imshow(Txyplt, interpolation='none', origin='lower',
                  extent=([min(loniv * 180.0 / pi),
                           max(loniv * 180.0 / pi),
                           min(lativ * 180.0 / pi),
                           max(lativ * 180.0 / pi)]))
plt.pyplot.xlabel('Longitude (deg)')
plt.pyplot.ylabel('Latitude (deg)')
plt.pyplot.title('Example No. 1: Synthesis from SRBFs\n'
                 '(Txy in E, spectral band %d -- %d)' % (nmin, nmax))
plt.pyplot.colorbar(orientation='horizontal')
plt.pyplot.grid()
plt.pyplot.show()
# -----------------------------------------------


# Read reference data for the element vxy that serve as sample reference
# data for the validation
# -----------------------------------------------
print('Loading reference data (the Txy element)')

pathtmp = os.path.join(path, 'Example1_ref_data.txt')
Txyref = np.loadtxt(pathtmp)
# -----------------------------------------------


# Perform the validation
# -----------------------------------------------
diffs = Txy * 1e9 - Txyref

print('')
print('Statistics of the differences between the SRBF synthesis and '
      'the reference values')

print('Min. (E) %0.16e' % min(diffs))
print('Max. (E) %0.16e' % max(diffs))
print('Mean (E) %0.16e' % np.mean(diffs))
print('STD (E) %0.16e' % np.std(diffs))

print('')
print('(These values should be generally below the 10**(-13) E level)')
# -----------------------------------------------


# Plot the validation results
# -----------------------------------------------
print('Plotting the validation results...')

diffs = diffs.reshape(nlati, nloni, order='F')
plt.pyplot.figure()
plt.pyplot.imshow(diffs, interpolation='none', origin='lower',
                  extent=([min(loniv * 180.0 / pi),
                           max(loniv * 180.0 / pi),
                           min(lativ * 180.0 / pi),
                           max(lativ * 180.0 / pi)]))
plt.pyplot.xlabel('Longitude (deg)')
plt.pyplot.ylabel('Latitude (deg)')
plt.pyplot.title('Example No. 1: SRBFs synthesis minus reference values\n'
                 '(Txy in E, spectral band %d -- %d)' % (nmin, nmax))
plt.pyplot.colorbar(orientation='horizontal')
plt.pyplot.grid()
plt.pyplot.show()
# -----------------------------------------------

print('===================================\n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# End of TEST EXAMPLE NO. 1





# TEST EXAMPLE NO. 2
# In this example, the expansion coefficients from "Test Example No. 1" are
# used to synthesize the Tzz element of the disturbing tensor, but this time
# within a difference spectral band than that for which the coefficients
# where estimated. More specifically, while the coefficients cover the
# harmonic degrees "nmin" -- "nmax", here, we use them to synthesize Tzz
# within the spectral band "nmin1" -- "nmax1", where "nmin < nmin1 <= nmax"
# and "nmin1 <= nmax1 <= nmax". The Tzz values synthesized from SRBFs are
# then compared with SHS as a reference. The validation is conducted on a
# sphere being 1000 m above the reference sphere.
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
print('')
print('')
print('===================================')
print('TEST EXAMPLE NO. 2\n')

nmin1 = 10  # Minimum degree of the synthesis for EXAMPLE NO. 2
nmax1 = 20  # Maximum degree of the synthesis for EXAMPLE NO. 2


# Perform the SRBFs synthesis
# -----------------------------------------------
print('SRBFs synthesis...')

(T, Tx, Ty, Tz, Txx, Txy, Txz, Tyy, Tyz, Tzz) = SP.srbfs_synthesis(lati,
 loni, ri,  lat, lon, r, nmin1, nmax1, a, phisrbf, R, psimax)
# Note that the syntax of "SP.srbfs_synthesis" is the same as in the MATLAB
# counterpart of the toolbox, but slightly different than in the Fortran
# subroutine "SRBFs_synthesis.f95". Here, three optional input parameters
# can be specified, defining the number of nodal points, the number of
# evaluation points and the size of the "phisrbf" vector.
#
# Using these optional parameter, the previous command reads
#
# (T, Tx, Ty, Tz, Txx, Txy, Txz, Tyy, Tyz, Tzz) = SP.srbfs_synthesis(lati,
# loni, ri,  lat, lon, r, nmin1, nmax1, a, phisrbf, R, psimax, npi, np2,
# len(phisrbf))
#
# and yields the same results.
# -----------------------------------------------


# Plot the SRBF synthesis
# -----------------------------------------------
print('Plotting the SRBF synthesis...')

Tzztmp = Tzz.reshape(nlati, nloni, order='F') * 1e9
plt.pyplot.figure()
plt.pyplot.imshow(Tzztmp, interpolation='none', origin='lower',
                  extent=([min(loniv * 180.0 / pi),
                           max(loniv * 180.0 / pi),
                           min(lativ * 180.0 / pi),
                           max(lativ * 180.0 / pi)]))
plt.pyplot.xlabel('Longitude (deg)')
plt.pyplot.ylabel('Latitude (deg)')
plt.pyplot.title('Example No. 2: SRBFs synthesis\n'
                 '(Tzz in E, spectral band %d -- %d)' % (nmin1, nmax1))
plt.pyplot.colorbar(orientation='horizontal')
plt.pyplot.grid()
plt.pyplot.show()
# -----------------------------------------------


# Read reference data for the element vxy that serve as sample reference
# data for the validation
# -----------------------------------------------
print('Loading reference data (the Txy element)')

pathtmp = os.path.join(path, 'Example2_ref_data.txt')
Tzzref = np.loadtxt(pathtmp)
# -----------------------------------------------


# Perform the validation
# -----------------------------------------------
diffs = Tzz * 1e9 - Tzzref

print('')
print('Statistics of the differences between the SRBF synthesis and '
      'the reference values')

print('Min. (E) %0.16e' % min(diffs))
print('Max. (E) %0.16e' % max(diffs))
print('Mean (E) %0.16e' % np.mean(diffs))
print('STD (E) %0.16e' % np.std(diffs))

print('')
print('(These values should be generally below the 10**(-12) E level)')
# -----------------------------------------------


# Plot the validation results
# -----------------------------------------------
print('Plotting the validation results...')

diffstmp = diffs.reshape(nlati, nloni, order='F')
plt.pyplot.figure()
plt.pyplot.imshow(diffstmp, interpolation='none', origin='lower',
                  extent=([min(loniv * 180.0 / pi),
                           max(loniv * 180.0 / pi),
                           min(lativ * 180.0 / pi),
                           max(lativ * 180.0 / pi)]))
plt.pyplot.xlabel('Longitude (deg)')
plt.pyplot.ylabel('Latitude (deg)')
plt.pyplot.title('Example No. 2: SRBFs synthesis minus reference values\n'
                 '(Tzz in E, spectral band %d -- %d)' % (nmin1, nmax1))
plt.pyplot.colorbar(orientation='horizontal')
plt.pyplot.grid()
plt.pyplot.show()
# -----------------------------------------------

print('===================================\n')
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# End of TEST EXAMPLE NO. 2



print('')
print('')
print('End of the test run.')
# ===============================================
