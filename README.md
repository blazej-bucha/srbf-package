# Introduction

`SRBF_package` is a small collection of a few Fortran and MATLAB routines for
gravity field modelling with band-limited spherical radial basis functions
(SRBFs), allowing both analysis and synthesis of the gravitational field.


## Features

* Harmonic analysis of any continuous scalar function with band-limited
  spherical radial basis functions.

* Synthesis of

  * the gravitational potential (or any other band-limited continuous function
    on a 2D sphere),

  * the gravitational vector in the local north-oriented reference frame
    (LNOF), and

  * the gravitational tensor in LNOF.

* Python wrapper available for the Fortran code.


# Documentation

The documentation can be found in
[docs/srbf-package.pdf](docs/srbf-package.pdf).



# Contact

Feel free to contact the author, Blazej Bucha, at blazej.bucha@stuba.sk.


# Citing

If you use `SRBF_package` in your work, please consider citing the associated
paper.

* Bucha, B., Janák, J., Papčo, J., Bezděk, A., 2016. High-resolution regional
  gravity field modelling in a mountainous area from terrestrial gravity
  data. Geophysical Journal International 207, 949-966,
  [http://doi.org/10.1093/gji/ggw311](http://doi.org/10.1093/gji/ggw311)


# Other related projects

* [CHarm](https://github.com/blazej-bucha/charm): C library for spherical
  harmonic transforms up to high degrees (tens of thousands and beyond).
  Supports OpenMP parallelization for shared memory architectures and
  vectorized CPU instructions (AVX, AVX2, AVX-512).

* [GrafLab](https://github.com/blazej-bucha/graflab) (GRAvity Field
  LABoratory): GrafLab (GRAvity Field LABoratory) is a MATLAB-based routine to
  compute functionals of the geopotential up to high degrees (tens of thousands
  and beyond).

* [isGrafLab](https://github.com/blazej-bucha/isgraflab) (Irregular Surface
  GRAvity Field LABoratory): A modified version of GrafLab to perform accurate
  and fast synthesis of gravity field quantities at dense grids residing on
  irregular surfaces such as the Earth's surface.

