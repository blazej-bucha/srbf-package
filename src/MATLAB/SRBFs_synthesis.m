function [V, Vx, Vy, Vz, Vxx, Vxy, Vxz, Vyy, Vyz, Vzz] = SRBFs_synthesis(...
    lati, loni, ri, lat, lon, r, nmin, nmax, a, phi, R, psimax)
% =========================================================================
%
% DESCRIPTION: This function synthesizes the gravitational potential, the
%              gravitational vector and the gravitational tensor from a given
%              set of expansion coefficients related to band-limited spherical
%              radial basis functions (SRBFs). 
%
%              For the adopted definition of SRBFs, see the attached PDF 
%              "Gravitational_potential_up_to_its_second-order_derivatives_in_terms_of_SRBFs.pdf". 
%              The gravitational vector and tensor are expressed in the local
%              north-oriented reference frame (LNOF). The definition of LNOF
%              as well as of the formulae used for the synthesis are provided
%              in the attached PDF.
%
%              The SRBF synthesis is parallelized over evaluation points via
%              the "parfor" loop. To enable parallel computation, open a parpool
%              session in the Command Window (command: parpool).
%
% INPUTS: "lati"   -- Spherical latitudes in radians of nodal points, at which
%                     SRBFs are formally located; a vector of dimensions (I, 1)
%
%         "loni"   -- Spherical longitudes in radians of nodal points, at which
%                     SRBFs are formally located; a vector of dimensions (I, 1)
%
%         "ri"     -- Spherical radius (in metres) of the sphere, at which SRBFs are
%                     formally located; a scalar ("ri" is assumed to be constant for
%                     all nodal points)
%
%         "lat"    -- Spherical latitudes (in radians) of evaluation points;
%                     a vector of dimensions (J, 1)
%
%         "lon"    -- Spherical longitudes (in radians) of evaluation points;
%                     a vector of dimensions (J, 1)
%
%         "r"      -- Spherical radii in metres of evaluation points;
%                     a vector of dimensions (J, 1)
%
%         "nmin"   -- Minimum harmonic degree of the synthesis; a scalar
%                     satisfying the inequality 0 <= nmin <= nmax
%
%         "nmax"   -- Maximum harmonic degree of the synthesis; a non-negative
%                     scalar
%
%         "a"      -- Expansion coefficients of SRBFs in m^4 * s^-2 (in the
%                     case of the input gravitational potential in m^2 * s^-2);
%                     a vector of dimensions (I, 1)
%
%         "phi"    -- Shape coefficients of SRBFs (dimensionless); a vector of 
%                     dimensions (nmax+1, 1)
%
%         "R"      -- Radius of the reference sphere, at which the nodal
%                     points (ri, lati, loni) are placed (metres); a scalar
%
%         "psimax" -- Spherical distance in radians, up to which SRBFs
%                     placed at the nodal points (lati, loni, ri) are
%                     considered in the synthesis; beyond psimax, the
%                     expansion coefficients of SRBFs at the nodal points 
%                     (lati, loni, ri) are automatically assumed to be equal
%                     to zero; for global integration, choose "psimax = pi",
%                     for regional (cap) integration, choose "0 < psimax < pi";
%                     a scalar
%
% OUTPUTS: "V"   -- Gravitational potential synthesized at the evaluation
%                   points (lat, lon, r) in m^2 * s^-2; a vector of dimensions (J, 1)
%
%          "Vx, Vy, Vz" -- Elements of the gravitational vector in LNOF
%                   synthesized at the evaluation points (lat, lon, r) in
%                   m * s^-2; three vectors of dimensions (J, 1)
%
%          "Vxx, Vxy, Vxz, Vyy, Vyz, Vzz" -- Elements of the gravitational tensor
%                   in LNOF synthesized at the evaluation points (lat, lon, r)
%                   in s^-2; six vectors of dimensions (J, 1)
%
% REFERENCES: Bucha, B., Bezdek, A., Sebera, J., Janak, J., 2015. Global and
%                regional gravity field determination from GOCE kinematic orbit
%                by means of spherical radial basis functions. Surveys in 
%                Geophysics 36, 773-801, http://doi.org/10.1007/s10712-015-9344-0
%                (preprint freely available at https://www.researchgate.net/profile/Blazej_Bucha)
%
%             Bucha, B., Janak, J., Papco, J., Bezdek, A., 2016. High-resolution
%                regional gravity field modelling in a mountainous area 
%                from terrestrial gravity data. Geophysical Journal 
%                International 207, 949-966, http://doi.org/10.1093/gji/ggw311
%                (preprint freely available at https://www.researchgate.net/profile/Blazej_Bucha)
%
%             Bucha, B., Hirt, C., Kuhn, M., 2019. Cap integration in spectral 
%                gravity forward modelling up to the full gravity tensor. 
%                Journal of Geodesy, https://doi.org/10.1007/s00190-019-01277-3
%                (preprint freely available at https://www.researchgate.net/profile/Blazej_Bucha)
%
%             Sprlak, M., Hamackova, E., Novak, P., 2015 Alternative
%                validation method of satellite gradiometric data by integral
%                transform of satellite altimetry data. Journal of Geodesy 89,
%                757-773, http://doi.org/10.1007/s00190-015-0813-5
%
% Contact: blazej.bucha@stuba.sk
%
%
% Please use the following reference when using this function:
%
%    Bucha, B., Janak, J., Papco, J., Bezdek, A., 2016. High-resolution
%        regional gravity field modelling in a mountainous area 
%        from terrestrial gravity data. Geophysical Journal 
%        International 207, 949-966, http://doi.org/10.1093/gji/ggw311
%
%
% Code history: Version 1.0 (Aug 9, 2019)
%
%                           -- The first published version of the code
%
%               Version 1.1 (Feb 18, 2020)
%
%                           -- The input vector "phi" now requires
%                              dimensions of (nmax + 1, 1) instead of 
%                              (1, nmax + 1)
%
%                           -- Minor formal modifications (improved 
%                              description, etc.)
%
% =========================================================================

if (length(lat) ~= length(lon)) || (length(lat) ~= length(r))  || (length(lon) ~= length(r))
    error('The vectors defining the position of evaluation points ("lat", "lon", "r") must be of the same dimensions.')
end

if length(lati) ~= length(loni)
    error('The vectors defining the position of nodal points ("lati", "loni") must be of the same dimensions.')
end

if length(lati) ~= length(a)
    error('The vector defining the expansion coefficients ("a") must be of the same dimensions as "lati" and "loni".')
end

if length(phi) ~= (nmax + 1)
    error('The vector defining shape coefficients of SRBFs ("phi") must have the size of (nmax + 1, 1).')
end


J = length(lat);   % Total number of evaluation points
I = length(lati);  % Total number of nodal points


% Transformation of spherical coordinates to geocentric cartesian
% coordinates
% --------------------------------------------------------------------------
[Xi, Yi, Zi] = sph2cart(loni, lati); % Unit vectors of the nodal points
[X, Y, Z] = sph2cart(lon, lat);      % Unit vectors of the evaluation points
% --------------------------------------------------------------------------


% Initializations
% --------------------------------------------------------------------------
V = zeros(J, 1); Vx = V; Vy = V; Vz = V; Vxx = V; Vxy = V; Vxz = V; 
Vyy = V; Vyz = V; Vzz = V;
% --------------------------------------------------------------------------


% Some useful substitutions
% --------------------------------------------------------------------------
degrees = 0:nmax; % Vector of harmonic degrees
R2 = R^2;
R2pi4 = 4 * pi * R2;
a = a';
phi = phi';
cospsimax = cos(psimax);
% --------------------------------------------------------------------------

phi(1:nmin) = 0; % If nmin>0, the shape coefficients of SRBFs are set to zero
                 % for correct SRBF synthesis

if abs(pi - psimax) < 1e-14 % Global integration domain
    int_domain = 0;
else % Regional integration domain
    int_domain = 1;
end

parfor j = 1:J % Loop over evaluation points
    
    r_j = r(j);
    lat_j = lat(j);
    lon_j = lon(j);
    
    cospsi_j = X(j) .* Xi + Y(j) .* Yi + Z(j) .* Zi; % Cosine of spherical
    % distances between the jth evaluation point and the nodal points
    
    if int_domain == 1 % Regional integration domain
        % Delete from the synthesis nodal points
        % located beyond the spherical distance "psimax". In other words, the
        % expansion coefficients for SRBFs placed at these nodal points are set
        % to zero; their contribution is therefore zero.
        idx = (cospsi_j >= cospsimax);
        cospsi_j = cospsi_j(idx); 
        lat_i_temp = lati(idx); %#ok<PFBNS>
        lon_i_temp = loni(idx); %#ok<PFBNS>
        I_idx = sum(idx);
        a_idx = a(idx); %#ok<PFBNS>
    else
        lat_i_temp = lati;
        lon_i_temp = loni;
        I_idx = I;
        a_idx = a;
    end
        
    cospsi_j(cospsi_j > 1) = 1; % To avoid numerical inaccuracies causing 
                                % imaginary numbers in the output quantities
    cospsi_j(cospsi_j < -1) = -1; % To avoid numerical inaccuracies causing 
                                  % imaginary numbers in the output quantities
    sinpsi_j = sin(acos(cospsi_j));
    sinpsi_j2 = sinpsi_j .* sinpsi_j;
    
    % Azimuth between nodal points and the jth evaluation point
    azimuth_j = atan2(cos(lat_i_temp) .* sin(lon_i_temp - lon_j), ...
        (cos(lat_j) .* sin(lat_i_temp) - sin(lat_j) .* cos(lat_i_temp).* ...
        cos(lon_i_temp - lon_j)));
    cosazimuth_j = cos(azimuth_j);
    sinazimuth_j = sin(azimuth_j);
    cos2azimuth_j = cos(2 * azimuth_j);
    sin2azimuth_j = sin(2 * azimuth_j);
    
    
    % Substitutions for Phi
    % ----------------------------------------------------------------------
    t = ri / r_j;
    cn = (2 * degrees + 1) ./ R2pi4;
    
    Phi_n = cn .* t.^(degrees+1) .* phi;
    % ----------------------------------------------------------------------
    
    
    % Substitutions for Phi10 and Phi11
    % ----------------------------------------------------------------------    
    Phi11_n = -1 ./ R .* t .* Phi_n;
    Phi10_n = -1 ./ R .* (degrees + 1) .* t .* Phi_n;
    % ----------------------------------------------------------------------
    
    
    % Substitutions for  Phi20, Phi21 and Phi22
    % ----------------------------------------------------------------------   
    tt = t .* t;
    
    Phi22_n =  1 ./ (2 * R2) .* tt .* Phi_n;
    Phi21_n = -1 ./ R2 .* (degrees+2) .* tt .* Phi_n;
    Phi20_n =  1 ./ R2 .* (degrees+1) .* (degrees + 2) .* tt .* Phi_n;
    % ----------------------------------------------------------------------
    
    
    % Computation of SRBFs Phi, Phi10, Phi11, Phi20, Phi21, Phi22
    % ----------------------------------------------------------------------
    % Seed values for Legendre polynomials (LPs) and their derivatives (dLP, ddLP)
    % ......................................................................
    P0 = ones(I_idx, 1);
    P1 = cospsi_j;
    dP1 = P0;
    ddP1 = 0;
    % ......................................................................
    
    % Sum over harmonic degree 0 and 1
    % ......................................................................
    % Degree 0
    Phi   = Phi_n(1) * P0;
    Phi10 = Phi10_n(1) * P0;
    Phi20 = Phi20_n(1) * P0;
    
    % Degree 0 + degree 1
    if nmax>0
        Phi   = Phi + Phi_n(2) * P1;
        Phi10 = Phi10 + Phi10_n(2) * P1;
        Phi11 = Phi11_n(2) * dP1;
        Phi20 = Phi20 + Phi20_n(2) * P1;
        Phi21 = Phi21_n(2) * dP1;
        Phi22 = 0;
    else
        Phi11 = 0;
        Phi21 = 0;
        Phi22 = 0;
    end
    % .......................;...............................................
    
    % Sum over harmonic degrees 2,...,nmax
    % ......................................................................    
    for n = 2:nmax
        % Recurrence relation for LPs
        P2 = ((2 * n - 1) / n) * cospsi_j .* P1 - ((n - 1) / n) * P0;
        
        % Recurrence relation for dLPs
        dP2 = n * P1 + cospsi_j .* dP1;
        
        % Recurrence relation for ddLPs
        ddP2 = (n + 1) * dP1 + cospsi_j .* ddP1;
        
        Phi   = Phi + Phi_n(n + 1) * P2;       
        Phi10 = Phi10 + Phi10_n(n + 1) * P2;
        Phi11 = Phi11 + Phi11_n(n + 1) * dP2;
        Phi20 = Phi20 + Phi20_n(n + 1) * P2;
        Phi21 = Phi21 + Phi21_n(n + 1) * dP2;
        Phi22 = Phi22 + Phi22_n(n + 1) * ddP2;
        
        P0 = P1;
        P1 = P2;
        dP1 = dP2;
        ddP1 = ddP2;
    end
    % ......................................................................
    Phi11 = Phi11 .* sinpsi_j;
    Phi21 = Phi21 .* sinpsi_j;
    Phi22 = Phi22 .* sinpsi_j2;
    % ----------------------------------------------------------------------

    % SRBF synthesis
    % ----------------------------------------------------------------------
    Vxx1 = a_idx * (-0.5 * Phi20);
    Vxx2 = a_idx * (cos2azimuth_j .* Phi22);
    
    V(j)   =  a_idx * Phi;
    Vx(j)  = -a_idx * (cosazimuth_j .* Phi11);
    Vy(j)  =  a_idx * (sinazimuth_j .* Phi11);
    Vz(j)  =  a_idx * Phi10;
    Vxx(j) =  Vxx1 + Vxx2;
    Vxy(j) = -a_idx * (sin2azimuth_j .* Phi22);
    Vxz(j) =  a_idx * (cosazimuth_j .* Phi21);
    Vyy(j) =  Vxx1 - Vxx2;
    Vyz(j) = -a_idx * (sinazimuth_j .* Phi21);
    Vzz(j) =  a_idx * Phi20;
    % ----------------------------------------------------------------------
end %End of the loop over evaluation points

end %End of the "SRBFs_synth" function
%==========================================================================


% SUBFUNCTIONS
% ==========================================================================
% --------------------------------------------------------------------------
function [X, Y, Z] = sph2cart(lon, lat)

% This function transforms spherical coordinates into cartesian coordinates,
% assuming that the point resides on the unit sphere.
%
% INPUT: "lon" -- vector of spherical longitudes in radians
%        "lat" -- vector of spherical latitudes in radians
%
% OUTPUT: "X","Y","Z" -- three vectors of cartesian coordinates

X = cos(lat) .* cos(lon);
Y = cos(lat) .* sin(lon);
Z = sin(lat);

end
% --------------------------------------------------------------------------
% ==========================================================================
