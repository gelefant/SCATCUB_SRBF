
function [XW,tri_vertices,tri_conn_list]=cub_sphpgon(n,vertices)

%--------------------------------------------------------------------------
% Input:
% n        : algebraic degree of exactness of the rule;
% vertices : it is a M x 3 matrix, where the k-th row represents the
%            cartesian coordinates of the k-th vertex of the spherical
%            polygon (counterclockwise order).
%--------------------------------------------------------------------------
% Output:
% XW       : N x 4 matrix, where the XW(k,1:3) are the cartesian
%            coordinates of the k-th node, while XW(k,4) is the respective
%            weight.
% tri_vertices,tri_conn_list: triangulation properties.
%--------------------------------------------------------------------------
% Subroutines:
% 1. triangulate_sphpgon_tg (attached);
% 2. cub_sphptri (external) and its subroutines.
%--------------------------------------------------------------------------
% Information:
% Authors: A. Sommariva and M. Vianello.
% Date: June 13, 2021.
%--------------------------------------------------------------------------
%% Copyright (C) 2021
%% Alvise Sommariva, Marco Vianello.
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
%% Alvise Sommariva, Marco Vianello.
%%
%% Date: JULY 13, 2021
%--------------------------------------------------------------------------


% ................. Troubleshooting and settings ..........................

if norm(vertices(1,:)-vertices(end,:))==0,vertices=vertices(1:end-1,:);end

% ........ Procedure start up, data stored in structures ............

%..........................................................................
% Main strategy:
% First we triangulate the spherical polygon, providing the vertices of the
% triangulation and its connectivity list (description of each spherical
% triangle via its vertices).
% Next we determine a rule for each spherical triangle.
%..........................................................................


%[tri_vertices,tri_conn_list]=triangulate_sphpgon(vertices);

subdom_index=find(isnan(vertices(:,1)) == 1);

% ....... Regularize if some vertex coordinates are NaNs ....... 

%..........................................................................
% NaN components are usually introduced in vertices components to separate 
% disconnected domains.
% We take this into account separating the analysis of the k-th subdomain
% taking data from vertices index "first_index(k)" to "last_index(k)".
%..........................................................................
first_index=1;
if length(subdom_index) == 0
    first_index=1;
    last_index=size(vertices,1);
else
    if subdom_index(end) == size(vertices,1)
        subdom_index=subdom_index(1:end-1);
    end
    
    first_index=[1 subdom_index+1];
    last_index=[subdom_index-1 size(vertices,1)];
end



% Computation of cubature rule over spherical polygon.
XW=[];
for k=1:length(first_index)
    
    vertices_sub=vertices(first_index:last_index,:); 
    [tri_vertices,tri_conn_list]=triangulate_sphpgon_tg(vertices_sub);
   
    for k=1:size(tri_conn_list,1)
        verticesL=tri_vertices((tri_conn_list(k,:))',:);
        P1=verticesL(1,:); P2=verticesL(2,:); P3=verticesL(3,:);
        XWL = cub_sphtri(n,P1,P2,P3); XW=[XW; XWL];
%         size(XWL,1)
    end
    
end








function [tri_vertices,tri_conn_list]=triangulate_sphpgon_tg(vertices)

%--------------------------------------------------------------------------
% Object:
% In this routine we determine a triangulation of the spherical polygon
% with M vertices described in cartesian coordinates the M x 3 matrix
% "vertices".
% The points of the triangulation are stored in "tri_vertices", while the
% K x 3 connectivity matrix is described in "tri_conn_list".
% In particular the k-th row of "tri_conn_list" consists of the indices of
% vertices of the k-th spherical triangle, w.r.t. "tri_vertices".
%--------------------------------------------------------------------------

if size(vertices,1) == 3
    tri_vertices=vertices; tri_conn_list=[1 2 3]; return;
end

% .................... compute barycenter vertices ........................

R=norm(vertices(1,:)); vertices=vertices/R;
CC=mean(vertices); CC=CC/norm(CC);

% ................ rotation matrix centroid to north pole .................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1];
rotmat=R1*R2; inv_rotmat=rotmat';

% ........................ rotate vertices ................................

vertices_NP=(rotmat*vertices')';

% ....... map radially vertices to plane tangent to north pole ............

XX_NP=vertices_NP(:,1); YY_NP=vertices_NP(:,2); ZZ_NP=vertices_NP(:,3);
rat=1./ZZ_NP;

XX_NPm=rat.*XX_NP; YY_NPm=rat.*YY_NP; ZZ_NPm=ones(size(XX_NPm));

% ............ polyshape of projected polygon to North Pole ...............

PG = polyshape(XX_NPm,YY_NPm);

% ............ triangulation of projected polygon to North Pole ...........

tri = triangulation(PG);

tri_vertices_NPm=PG.Vertices;
tri_vertices_NPm=[tri_vertices_NPm ones(size(tri_vertices_NPm,1),1)];
rads=sqrt((tri_vertices_NPm(:,1)).^2+(tri_vertices_NPm(:,2)).^2+1);

% ......................... output data ...................................

% Once we have a triangulation at the North-Pole, the points with the same
% triangulation indices on the spherical polygon determine its
% triangulation.
tri_vertices_NP=tri_vertices_NPm./rads;
tri_vertices=R*(inv_rotmat*tri_vertices_NP')';
tri_conn_list=tri.ConnectivityList;












