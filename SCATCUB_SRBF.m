function I = SCATCUB_SRBF(centers, fdata, domain, rbf_type, isdeg,verbose)
% -------------------------------------------------------------------------
% Function for the algorithm described in the paper 
% R.Cavoretto, A.De Rossi, G.Elefante and A.Sommariva 
% "Adaptive RBF cubature by scattered data on spherical polygons"
%
% Starting from scattered data, the integral of the underlying function is 
% computed by approximating it with an RP or TPS RBF interpolant constructed 
% on the original data. The exponent parameter is chosen through a variant 
% of LOOCV that allows selecting the best one for integration.
% -------------------------------------------------------------------------
% INPUT 
% centers  : scattered points in the spherical polygon
% fdata    : values of the underlying function in the scattered points
% domain   : polyshape of the spherical polygon
% rbf_type : setting the rbf - 1: RP; 2: TPS
% isdeg    : expressing if the vertices of the domain are in deg or rad
% verbose  : 1: text; 2: no text
% -------------------------------------------------------------------------
% OUTPUT
% I        : approximation of the integral via the rbf chosen and an
%            adaptive integration
%--------------------------------------------------------------------------
%% Copyright (C) 2021-
%% Roberto Cavoretto, Alessandra De Rossi, Giacomo Elefante, Alvise Sommariva.
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
%%
%% Authors: 
%% Roberto Cavoretto, Alessandra De Rossi, Giacomo Elefante, Alvise Sommariva.
%--------------------------------------------------------------------------

if nargin < 4 || isempty(rbf_type)
    rbf_type = 2;
end

if nargin < 5
    isdeg = 1; % 1 if coordinates are in degree, 0 if are in radiant
end

if nargin < 6
    verbose = 1;
end

% Extracting the vertices of the domain
vertices = domain.Vertices;

if isdeg == 1
    % Mapping coordinates in radiant
    [Vx,Vy,Vz] = sph2cart(deg2rad(vertices(:,1)),deg2rad(vertices(:,2)),1);
    vertices = [Vx,Vy,Vz];
end

switch rbf_type
    case 1
        phi = @(n,r) (2-2*r).^((2*n-1)/2);   % Radial power (RP)
        % Parameters for RBF (power) parameter loop below
        minn = 1; maxn = 8;
        n = minn:maxn;
    case 2
        phi = @(n,r) (1-r).^(n).*log(2-2*r+10^(-50)); % Thin-Plate Spline (TPS)
        % Parameters for RBF (power) parameter loop below
        minn = 1; maxn = 8;
        n = minn:maxn;
end

% Computing the area of the spherical polygon
XWC=cub_sphpgon(1,vertices);
area = sum(XWC(:,4));

DM_data = centers*centers';

% Choice of the best exponent via a modified-LOOCV for integration
avgEF = Inf*ones(size(n));
condV = zeros(size(n));
for i = 1:length(n)
    IM=phi(n(i),DM_data);
    condV(i) = condest(IM);
    if condV(i)>1/eps %1/eps %1e+20 
        continue
    end
    invIM = pinv(IM);
    coeff=IM\fdata;

    EF = (coeff)./diag(invIM);

    % Compute the average of the norm of EF
    avgEF(i) = abs(sum(area/size(centers,1)*EF(:)));
end
[minAvgEF,idn] = min(avgEF);
% save V1000.mat condV
if verbose
    fprintf('\n \n \t RBF EXP CHOSEN             : %d ',n(idn));
end

% Computing the coefficients of the RBF-interpolant
V=phi(n(idn),centers*centers');
coeff=V\fdata; 

% Computing the integration of the RBF-interpolant
atol=10^(-12);
rtol=atol;
I = adaptive_RBF_sphpgon(vertices,coeff,phi,n(idn),centers,atol,rtol);

end
