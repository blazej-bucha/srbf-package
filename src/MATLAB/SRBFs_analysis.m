function [a, lati, loni] = SRBFs_analysis(lati, loni, T, R, nmax)
% =========================================================================
%
% DESCRIPTION: This function performs surface SRBF analysis up to degree "nmax"
%              of a function sampled at the Driscoll--Healy grid residing 
%              on the reference sphere with the radius "R".
%
%
% INPUTS: "lati"   -- Vector of spherical latitudes (in radians) over a single 
%                     meridian of the Driscoll--Healy grid (the first output from the
%                     "DH_grid" function); dimensions of the vector (2 * (nmax + 1) + 1, 1)
%
%         "loni"   -- Vector of spherical longitudes (in radians) over a single 
%                     parallel of the Driscoll--Healy grid (the second output from the
%                     "DH_grid" function); dimensions of the vector (2 * (nmax + 1), 1)
%
%         "T"      -- Input data that are analysed; matrix of dimensions
%                     (2 * (nmax + 1) + 1, 2 * (nmax + 1)) with the following
%                     structure
%
%                     [T(lati(1),loni(1))            T(lati(1),loni(2))            ...            T(lati(1),loni(2*(nmax+1),1))]
%                     [T(lati(2),loni(1))            T(lati(2),loni(2))            ...            T(lati(2),loni(2*(nmax+1),1))]
%                     [       .                                                    .                         .                 ]
%                     [       .                                                     .                        .                 ]
%                     [       .                                                      .                       .                 ]            
%                     [T(lati(2*(nmax+1)+1),loni(1)) T(lati(2*(nmax+1)+1),loni(2)) ... T(lati(2*(nmax+1)+1),loni(2*(nmax+1),1))]
%
%                     where "lati" and "loni" denote the output vectors from
%                     the "DH_grid" function.
%
%         "R"      -- Radius of the reference sphere (metres) at which the data are
%                     sampled and SRBFs are formally located; the output expansion
%                     coefficients "a" refer to this sphere; a scalar
%
%         "nmax"   -- Maximum harmonic degree up to which the analysis is
%                     performed, nmax > 0; a scalar
%
%
% OUTPUTS: "a"    -- Expansion coefficients; units are equal to the units
%                    of the input data but multiplied by m^2 in addition)
%                    at the nodal points (lati,loni,ri); dimensions of the vector 
%                    ((2 * (nmax + 1) + 1) * (2 * (nmax + 1)), 1)
%
%          "lati" -- Spherical latitudes in radians of the nodal points, to
%                    which expansion coefficients refer to; dimensions of 
%                    the vector ((2 * (nmax + 1) + 1) * (2 * (nmax + 1)), 1)
%
%          "loni" -- Spherical longitudes in radians of the nodal points, to
%                    which expansion coefficients refer to; dimensions of 
%                    the vector (( 2 * (nmax + 1) + 1) * (2 * (nmax + 1)), 1)
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
colati = pi / 2 - lati;

% Integration weights
% --------------------------------------------------------------------------
w = 0;
for k = 0:(L - 1)
    w = w + 1 / (2 * k + 1) * sin((2 * k + 1) *colati);
end
w = 2 * pi * R^2 / L^2 * sin(colati) .* w;
%--------------------------------------------------------------------------


% SRBF analysis
a = bsxfun(@times, w, T);
a = a(:);


[loni, lati] = meshgrid(loni, lati); % Generates latitudes and longitudes of
% all nodal points, so that the vectors "lati" and "loni" match the
% expansion coefficients from the vector "a"
loni = loni(:);
lati = lati(:);

end
