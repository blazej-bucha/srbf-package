% This is a script to test the SRBF_package. The script consists of
% two examples that demonstrate the correctness of the package
% in a closed-loop environment.
% The SRBF synthesis is parallelized over evaluation points via
% the "parfor" loop. To enable parallel computation, open a parpool
% session in the Command Window (command: parpool).

clear;
clc;
close all;


% Inputs
% ==========================================================================
nmin = 0;  % Minimum harmonic degree of the expansion
nmax = 35; % Maximum harmonic degree of the expansion

R = 6378136.3; % Radius of the reference sphere (has to be equal to the radius of the
             % nodal points "ri"; see below)


% Nodal points
% --------------------------------------------------------------------------
[lati, loni] = DH_grid(nmax); % Generates latitudes and longitudes (in radians)
                           % of the Driscoll--Healy grid points for current "nmax"
ri = R; % Spherical radius (in metres) of the nodal points (has to be the same
      % for all nodal points)
%--------------------------------------------------------------------------


% Evaluation points
% --------------------------------------------------------------------------
% In this example, the evaluation points are the same as the nodal points,
% but this may be changed. The radius of evaluation points is 1000 m above
% the reference sphere, at which the nodal points are placed.
[lat, lon] = DH_grid(nmax); % Latitudes and longitudes of the evaluation points
                            % (in radians)
r = (R + 1000); % Spherical radius of the evaltuation points
                % (may be different for each point, but the inequality "r>=R" 
                % must be satisfied)
%--------------------------------------------------------------------------

psimax = pi; % Integration radius in radians (for global integration, choose
             % "psimax = pi"; for regional (cap) integration, choose
             % 0 < psimax < pi)
phi = ones(nmax + 1, 1); % Shape coefficients of SRBFs. In this example, we use
                         % the Shannon SRBFs with shape coefficients being equal
                         % to one up to the "nmax" degree and zero beyond.
%==========================================================================


nlati = length(lati); % Number of latitudes at one meridian of the nodal points grid
nloni = length(loni); % Number of longitudes at one parallel of the nodal points

nlat = length(lat); % Number of latitudes at one meridian of the grid with evaluation points
nlon = length(lon); % Number of longitudes at one meridian of the grid with evaluation points

[lon_grd, lat_grd] = meshgrid(lon, lat); % Generates lats and lons for each evaluation point of the grid
lat_grd = lat_grd(:); % Transforms to column vector
lon_grd = lon_grd(:); % Transforms to column vector
r_grd = r * ones(length(lat_grd), 1); % Spherical radius of the evaluation points


% TEST EXAMPLE NO. 1
% In this example, the input quantity (the disturbing potential) is
% synthesized from global geopotential model EGM96 via spherical harmonic
% synthesis (SHS). The minimum and the maximum degrees of the synthesis are
% "nmin" and "nmax", respectively (specified in the "Inputs" section). The
% disturbing potential is synthesized on the reference sphere with the
% radius "R". These are the input data for SRBFs analysis.
% Through the SRBFs analysis, the expansion coefficients "a" are estimated via
% the exact quadrature.
% Finally, the Txy element of the disturbing tensor is synthesized
% from the expansion coefficients (via SRBF synthesis). To check the correctness
% of the upward continuation operations, the validation is performed on a 
% sphere being 1000 m above the reference sphere. These results are
% then compared with reference values from SHS. The differences between the
% two results reveal the associated errors of the SRBFs analysis and
% synthesis (errors in SHS are here neglected).
% ==========================================================================
fprintf('TEST EXAMPLE NO. 1\n')
fprintf('================================\n')
fprintf('Computing input quantity for SRBFs analysis... (%s)\n', datestr(clock))
%--------------------------------------------------------------------------
T = GrafLab('OK', 3986004.415E+8, 6378136.3, nmin, nmax, 1, 'EGM96.mat', 1, ...
    0,lati * 180 / pi, 'empty', 'empty', loni * 180 / pi, 'empty', ...
    'empty', ri - 6378136.3, [], [], [], [], 'Out_file', 0, 5, 1, [], 0, ...
    0, 0, [], [], [], [], [], 0);
T = reshape(T(:, end), length(lati), length(loni));
%--------------------------------------------------------------------------


fprintf('Plotting input quantity for SRBFs analysis... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
figure;
imagesc(loni * 180 / pi, lati * 180 / pi, flipud(T));
colorbar
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title(sprintf('Example No. 1: Input quantity for the SRBFs analysis\n(the disturbing potential in m^2 s^{-2}, spectral band %d -- %d)', nmin, nmax))
caxis([-2.5 * std(T(:)), 2.5 * std(T(:))])
% --------------------------------------------------------------------------


fprintf('Computing reference values... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
Txyref = GrafLab('OK', 3986004.415E+8, 6378136.3, nmin, nmax, 1, 'EGM96.mat', ...
    1, 0, lat * 180 / pi, 'empty', 'empty', lon * 180 / pi, 'empty', ...
    'empty', r - 6378136.3, [], [], [], [], 'Out_file', 0, 9, 1, [], 0, ...
    0, 0, [], [], [], [], [], 0);
Txyref = reshape(Txyref(:, end-2), nlat, nlon);
% --------------------------------------------------------------------------


% SRBFs analysis -- Computation of the expansion coefficents via the numerical 
% quadrature (an exact method)
fprintf('SRBFs analysis... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
[a, lati_grd, loni_grd] = SRBFs_analysis(lati, loni, T, R, nmax);
% --------------------------------------------------------------------------


% Plot expansion coefficients
% --------------------------------------------------------------------------
fprintf('Plotting expansion coefficients... (%s)\n', datestr(clock))

figure;
scatter(loni_grd * 180 / pi, lati_grd * 180 / pi, 15, a, 'filled');
colorbar
xlabel('Longitide (deg)')
ylabel('Latitude (deg)')
title('Example No. 1: Expansion coefficients of SRBFs (m^4 s^{-2})')
xlim([min(loni) max(loni)] * 180 / pi);
ylim([min(lati) max(lati)] * 180 / pi);
box on
caxis([-2.5 * std(a(:)), 2.5 * std(a(:))])
% --------------------------------------------------------------------------


% SRBFs synthesis
fprintf('SRBFs synthesis... (%s)\n', datestr(clock))
%--------------------------------------------------------------------------
[T, Tx, Ty, Tz, Txx, Txy, Txz, Tyy, Tyz, Tzz] = SRBFs_synthesis(lati_grd, ...
    loni_grd, ri, lat_grd, lon_grd, r_grd, nmin, nmax, a, phi, R, psimax); %#ok<ASGLU>
% Note that 10 quantities are synthesized here, but only Txy is used in 
% Example 1.
Txy = reshape(Txy, nlati, nloni) * 10^9; %Transforms the column vector to a matrix
% and changes units to Eotvos (1 E = 10^-9 s^-2).
% --------------------------------------------------------------------------


fprintf('Plotting SRBFs synthesis... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
figure;
imagesc(loni * 180 / pi, lati * 180 / pi, flipud(Txy));
colorbar
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title(sprintf('Example No. 1: Synthesis from SRBFs\n(T_{xy} in E, spectral band %d -- %d)', nmin, nmax))
caxis([-2.5 * std(Txy(:)), 2.5 * std(Txy(:))])
% --------------------------------------------------------------------------


fprintf('Plotting differences between synthesized and reference values... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
figure;
imagesc(loni * 180/pi, lati * 180 / pi, flipud(Txy - Txyref));
colorbar
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title(sprintf('Example No. 1: SRBFs synthesis minus reference values\n(T_{xy} in E, spectral band %d -- %d)', nmin, nmax))
caxis([-2.5 * std(Txy(:) - Txyref(:)), 2.5 * std(Txy(:) - Txyref(:))])
% --------------------------------------------------------------------------
fprintf('================================\n\n')
%==========================================================================


% TEST EXAMPLE NO. 2
% In this example, the expansion coefficients from "Test Example No. 1" are
% used to synthesize the Tzz element of the disturbing tensor, but this time
% within a difference spectral band than that for which the coefficients where
% estimated. More specifically, while the coefficients cover the harmonic 
% degrees "nmin" -- "nmax", here, we use them to synthesize Tzz within 
% the spectral band "nmin1" -- "nmax1", where "nmin < nmin1 <= nmax" and 
% "nmin1 <= nmax1 <= nmax". The Tzz values synthesized from SRBFs are then
% compared with SHS as a reference. The validation is conducted on a sphere
% being 1000 m above the reference sphere.
% ==========================================================================
fprintf('TEST EXAMPLE NO. 2\n')
fprintf('================================\n')
nmin1 = 10; % Minumum degree of the synthesis
nmax1 = 20; % Maximum degree of the synthesis

phi = ones(nmax1+1, 1);

if (nmin1 <= nmin) || (nmin1 > nmax)
    error('Please choose an "nmin1" value that is larger than %d, but smaller than or equal to %d.', nmin, nmax)
end
if (nmax1 < nmin1) || (nmax1 > nmax)
    error('Please choose an "nmax1" value that is equal or larger than %d, but equal or smaller than %d.', nmin1, nmax)
end

% SRBFs synthesis
fprintf('SRBFs synthesis... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
[T, Tx, Ty, Tz, Txx, Txy, Txz, Tyy, Tyz, Tzz] = SRBFs_synthesis(lati_grd, ...
    loni_grd, ri, lat_grd, lon_grd, r_grd, nmin1, nmax1, a, phi, R, psimax);
% Note that 10 quantities are synthesized here, but only Txy is used in 
% Example 1.
Tzz = reshape(Tzz, nlati, nloni) * 10^9; %Transforms the column vector to a matrix
% and changes units to Eotvos (1 E = 10^-9 s^-2).
% --------------------------------------------------------------------------


fprintf('Plotting SRBFs synthesis... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
figure;
imagesc(loni * 180 / pi, lati * 180 / pi, flipud(Tzz));
colorbar
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title(sprintf('Example No. 2: SRBFs synthesis\n(T_{zz} in E, spectral band %d -- %d)', nmin1, nmax1))
caxis([-2.5 * std(Tzz(:)), 2.5 * std(Tzz(:))])
% --------------------------------------------------------------------------


% Reference values from SHS
fprintf('Computing reference values... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
Tzzref = GrafLab('OK', 3986004.415E+8, 6378136.3, nmin1, nmax1, 1, 'EGM96.mat', ...
    1, 0, lat * 180 / pi, 'empty', 'empty', lon * 180 / pi, 'empty', 'empty', ...
    r - 6378136.3, [], [], [], [], 'Out_file', 0, 8, 1, [], 0, 0, 0, [], [], ...
    [], [], [], 0);
Tzzref = reshape(Tzzref(:, end), length(lati), length(loni));
Tzzref = reshape(Tzzref, nlati, nloni);
% --------------------------------------------------------------------------


fprintf('Plotting differences... (%s)\n', datestr(clock))
% --------------------------------------------------------------------------
figure;
imagesc(loni * 180 / pi, lati * 180 / pi, flipud(Tzz - Tzzref));
colorbar
xlabel('Longitude (deg)')
ylabel('Latitude (deg)')
title(sprintf('Example No. 2: SRBFs synthesis minus reference values\n(T_{zz} in E, spectral band %d -- %d)', nmin1, nmax1))
caxis([-2.5 * std(Tzz(:) - Tzzref(:)), 2.5 * std(Tzz(:) - Tzzref(:))])
% --------------------------------------------------------------------------

fprintf('================================\n\n')

fprintf('End of the test run. (%s)\n', datestr(clock))
%==========================================================================
