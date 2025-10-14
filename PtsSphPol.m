function Y = PtsSphPol(N,pts_type,vertices)
% -------------------------------------------------------------------------
% This functions compute N points of three different type in a spherical 
% polygon of vertices "vertices".
% INPUT
% N         :  Points in the spherical cap
% pts_type  :  String for different sets of points, 
%              'H' are Halton points, 
%              'S' Sobol points,
%              'R' are random 
% vertices  :  vector Mx3 of the (cartesian) coordinates of the M vertices 
%              of the spherical polygon
% OUTPUT   
% Y         :  vector Nx3 of N points in the spherical polygon 
% -------------------------------------------------------------------------
% Dates
%--------------------------------------------------------------------------
% First version: June 15, 2024;
% Checked: November 27, 2024.
%--------------------------------------------------------------------------
%% Copyright (C) 2021-
%% Giacomo Elefante, Alvise Sommariva.
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
%% Giacomo Elefante, Alvise Sommariva.
%--------------------------------------------------------------------------
% .................... compute barycenter vertices ........................

R=norm(vertices(1,:)); vertices=vertices/R;
if  anynan(vertices(:,1))
    id = 1:size(vertices,1); idx = isnan(vertices(:,1)); nanidx = id(idx');
    CC=mean(vertices(1:(nanidx-1),:),1); CC=CC/norm(CC);
else
    CC=mean(vertices,1); CC=CC/norm(CC);
end

% ................ rotation matrix centroid to north pole .................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1];
rotmat=R1*R2; inv_rotmat=rotmat';

% ........................ rotate vertices to north pole...................

vertices_NP=(rotmat*vertices')';

% ..................... uniform points in the sphere cap...................

z_min = min(vertices_NP(:,3));

switch pts_type
    case 'H'
        p = haltonset(2);
    case 'S'
        p = sobolset(2);
end

% ................... stereographic map from south pole ...................

% ....... vertices

XX_SP=vertices_NP(:,1); YY_SP=vertices_NP(:,2); ZZ_SP=vertices_NP(:,3);
rat=1./(1+ZZ_SP);

XX_SPm=rat.*XX_SP; YY_SPm=rat.*YY_SP;

% ........... polyshape of projected polygon from South Pole ..............

PG = polyshape(XX_SPm,YY_SPm);

iter = 1;
Y = [];
while size(Y,1) < N
    Niter = N*iter;
    switch pts_type
        case 'R'
            Xnew = rand(Niter,2);
        otherwise
            Xnew = net(p,Niter);
    end
    Xnew = [z_min+(1-z_min)*Xnew(:,1),2*pi*Xnew(:,2)];
    Xnew3 = [sqrt(1-Xnew(:,1).^2).*cos(Xnew(:,2)), sqrt(1-Xnew(:,1).^2).*sin(Xnew(:,2)), Xnew(:,1)];

    rat2=1./(1+Xnew3(:,3));
    X_SP=rat2.*Xnew3(:,1); Y_SP=rat2.*Xnew3(:,2);
    in = isinterior(PG,X_SP,Y_SP);
    Xin = [X_SP(in),Y_SP(in)];
    if isempty(Xin)
        continue
    end
    XinS = [2*Xin(:,1)./(1+Xin(:,1).^2+Xin(:,2).^2),2*Xin(:,2)./(1+Xin(:,1).^2+Xin(:,2).^2),(1-Xin(:,1).^2-Xin(:,2).^2)./(1+Xin(:,1).^2+Xin(:,2).^2)];

    Y = (inv_rotmat*XinS')';
    iter = iter + 1;
end
Y = Y(1:N,:);

end