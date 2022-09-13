function [lati, loni] = DH_grid(nmax)
% =========================================================================
%
% DESCRIPTION: This function generates spherical latitudes and longitudes 
%              of the Driscoll--Healy grid points on the unit sphere (in radians).
%
%
% INPUT: "nmax" -- Maximum harmonic degree; a positive scalar
%
%
% OUTPUTS: "lati" -- Spherical latitudes (in radians) of the Driscoll--Healy 
%                    grid points at a single meridian; a vector of dimensions 
%                    (2 * (nmax + 1) + 1, 1)
%
%          "loni" -- Spherical longitudes (in radians) of the Driscoll--Healy
%                    grid points at a single parallel, a vector of dimensions 
%                    (2 * (nmax + 1), 1);
%
%
% REFERENCES: Driscoll, J. R., Healy, D.M., 1994. Computing Fourier transforms
%                and convolutions on the 2-sphere. Advances in Applied 
%                Mathematics 15, 202-250, doi: http://doi.org/10.1006/aama.1994.1008
%
%             Schmidt, M., Fengler, M., Mayer-Gurr, T., Eicker, A., Kusche,
%                J., 2007. Regional gravity modelling in terms of spherical
%                base functions. Journal of Geodesy 81, 17-38, 
%                doi: https://doi.org/10.1007/s00190-006-0101-5
%
%
% Contact: blazej.bucha@stuba.sk
%
%
% Code history: Version 1.0 (Aug 9, 2019)
%
%                           -- The first published version of the code
%
%               Version 1.1 (Feb 18, 2020)
%
%                           -- Minor formal modifications (improved 
%                              description, etc.)
%
% =========================================================================

L = nmax + 1;
j = 0:(2 * L);
k = 0:(2 * L - 1);

lati = -pi / 2 + pi .* j ./ (2 * L);
loni = pi * k ./ L;

lati = lati(:);
loni = loni(:);

end
