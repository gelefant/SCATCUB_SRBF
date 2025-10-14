
function [I,Ierr,flag,Ihigh,iters,tri_vertices,tri_conn_list,L1_vertices]=...
    adaptive_cub_sphpgon(vertices,f,atol,rtol)

%--------------------------------------------------------------------------
% Input:
% vertices : it is a M x 3 matrix, where the k-th row represents the
%            cartesian coordinates of the k-th vertex of the spherical
%            polygon (counterclockwise order);
% f        : function to integrate over the spherical triangle defined by
%            vertices;
% atol      : absolute cubature error tolerance.
% rtol      : relative cubature error tolerance.
%--------------------------------------------------------------------------
% Output:
% I      : approximation of integral of "f" over the spherical triangle
%           defined by vertices
%--------------------------------------------------------------------------
% Subroutines:
% 1. area_spherical_triangle;
% 2. center_spherical_triangle;
% 3. generate_sphtri_sons.
%--------------------------------------------------------------------------
% Information:
% Authors: A. Sommariva and M. Vianello.
% Date: June 7, 2021.
% Modified: November 12, 2024.
%--------------------------------------------------------------------------
%% Copyright (C) 2021-
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
%--------------------------------------------------------------------------

%fprintf('\n adaptive_cub_sphpgon');

multiple_domain = 0;

% ................. Troubleshooting and settings ..........................

if nargin < 3, atol=10^(-6); end
if nargin < 4, rtol=10^(-6); end

% max number of triangles in which the domain is partitioned
max_triangles=1000;

if norm(vertices(1,:)-vertices(end,:))==0,vertices=vertices(1:end-1,:);end

% ........ Procedure start up, data stored in structures ............

%..........................................................................
% Main strategy:
% First we triangulate the spherical polygon, providing the vertices of the
% triangulation and its connectivity list (description of each spherical
% triangle via its vertices).
%
% Below, we make 2 lists of sph. triangle data:
% 1. LIST 1: all the triangles provide contribution to the integration
%    result and for each triangle an error estimate is available;
% 2. LIST 2: all the triangles are "sons" triangles of some triangle in
%    LIST 1; integrals are computed in each of them but no error estimate
%    is available.
%
% We store the relevant data of each spherical triangle (LIST 1 set).
%
% * L1_vertices : cell, whose k-th element are the vertices of the k-th
%                 sph. triangle stored as a 3 x 3 matrix, whose j-th row
%                 are the cartesian coordinates  of the j-th vertex of the
%                 k-th sph. triangle (ordered counterclockwise);
% * L1_areas    : vector, whose k-component represent the area of the k-th
%                 sph. triangle;
% * L1_integrals: vector, whose k-component represent the approximation of
%                 the integral of "f" on the k-th sph. triangle;
% * L1_errors   : vector, whose k-component represent the approximation of
%                 the error of integral on the k-th sph. triangle.
% 
% If the error estimate is not satisfied, the code find the worst triangle,
% subdvides it, and updates the lists.
%..........................................................................

% fprintf('\n \t Adaptive code: triangulation');
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
elseif multiple_domain == 1
    if subdom_index(end) == size(vertices,1)
        subdom_index=subdom_index(1:end-1);
    end
    
    first_index=[1 subdom_index-1];
    last_index=[subdom_index+1 size(vertices,1)];
else
    first_index=1;
    last_index=size(vertices,1);
end



% Computation of cubature rule over spherical polygon.
XW=[];
tri_vertices=[]; tri_conn_list=[];
for k=1:length(first_index)
    
    vertices_sub=vertices(first_index:last_index,:); 
    [tri_verticesL,tri_conn_listL]=triangulate_sphpgon_tg(vertices_sub);
    tri_vertices=[tri_vertices; tri_verticesL]; 
    tri_conn_list=[tri_conn_list; tri_conn_listL];
end


L1_vertices={}; L1_integrals=[]; L1_errors=[]; L2_data={};

% fprintf('\n \t Adaptive code: start up');
for k=1:size(tri_conn_list,1)
    vertices=tri_vertices((tri_conn_list(k,:))',:);
    [L1_vertices,L1_integrals,L1_errors,L2_data,Itemp_low,Itemp_high]=...
        start_up(vertices,f,L1_vertices,L1_integrals,L1_errors,L2_data);
end

I=sum(L1_integrals);

I0=I;


% ........................ successive refinement ..........................

iters=1;


% fprintf('\n \t Adaptive code: refining');
while  sum(L1_errors) > max([atol rtol*I])
    
    % fprintf('\n \t %6.0f',iters);
    
    N=length(L1_integrals);
    
    % too many triangles: exit with errors
    if N > max_triangles
        I=sum(L1_integrals); Ierr=sum(L1_errors); flag=1;
        if nargout > 3, Ihigh=Ihigh_computation(L2_data); end
        return;
    end
    
    [Ierr_max,kmax]=max(L1_errors);
    
    % Erase "kmax" sph. triangle from LIST 1
    k_ok=setdiff(1:N,kmax);
    
    if length(k_ok) > 0
        L1_vertices={L1_vertices{k_ok}}; L1_integrals=L1_integrals(k_ok);
        L1_errors=L1_errors(k_ok);
    else
        L1_vertices={}; L1_integrals=[]; L1_errors=[];
    end
    
    % Move sub sph. triangles relative to "k-max" from LIST 2 to LIST 1.
    
    L2_data_kmax=L2_data{kmax}; % cell with 4 data structs
    
    if length(k_ok) > 0
        % erase cell relative to sub triangles of kmax
        L2_data={L2_data{k_ok}};
    else
        L2_data={};
    end
    % triangle, from LIST 2
    
    %......................................................................
    % Generate sons of each son of L2_data_kmax and compute approximate
    % integration errors for each son of L2_data_kmax.
    % 1. Each son of L2_data_kmax is moved to LIST 1.
    % 2. The son of each son of L2_data_kmax is moved to LIST 2.
    %......................................................................
    
    for j=1:4
        L2_data_temp=L2_data_kmax{j}; % struct data
        L2_data_temp_sons=generate_sphtri_sons(L2_data_temp.vertices,f);
        
        % "L2_data_kmax_j_sons" is a cell with 4 structs data.
        for jj=1:4
            L2_data_temp_sons_integral(jj)=L2_data_temp_sons{jj}.integral;
        end
        Itemp_low=L2_data_temp.integral;
        Itemp_high=sum(L2_data_temp_sons_integral);
        Itemp_err=abs(Itemp_high-Itemp_low);
        
        % update LIST 1 with j-th sub triangle "generated" by kmax.
        L1_vertices{end+1}=L2_data_temp.vertices;
        L1_integrals(end+1)=L2_data_temp.integral;
        L1_errors(end+1)=Itemp_err;
        
        % update LIST 2 with sons of j-th sub triangle "generated" by kmax.
        L2_data{end+1}=L2_data_temp_sons;
    end
    
    iters=iters+1;
    I=sum(L1_integrals);
    
end

% ....................... providing results ...............................

I=sum(L1_integrals); Ierr=sum(L1_errors); flag=0;
if nargout > 3, Ihigh=Ihigh_computation(L2_data); end

% fprintf('\n \t I0: %1.15e',I0);
% fprintf('\n \t I1: %1.15e',I);
















function [L1_vertices,L1_integrals,L1_errors,L2_data,Itemp_low,...
    Itemp_high]=start_up(vertices,f,L1_vertices,L1_integrals,L1_errors,...
    L2_data)

%--------------------------------------------------------------------------
% Object:
% Given a sph.triangle with vertices "vertices", a cell with previous
% vertices of previous triangles, a vector "L1_integrals" with integrals of
% "f" in the previous triangles and a vector "L1_errors" of the error
% estimates of the integrals of "f" in the previous triangles, it updates
% these cells and vectors with vertices of the new triangle, the integral
% of "f" in the new triangle, and its error estimate.
% Data about the sons of the new triangle is put in the cell "L2_data".
%--------------------------------------------------------------------------


% vertices
L1_vertices{end+1}=vertices;

% integrals
P1=(vertices(1,:))'; P2=(vertices(2,:))'; P3=(vertices(3,:))';

xyzw= cub_sphtri(4,P1,P2,P3);
Itemp_low=(xyzw(:,4))'*feval(f,xyzw(:,1),xyzw(:,2),xyzw(:,3));

% approximation of integral on triangle (few evals)
L1_integrals=[L1_integrals; Itemp_low];

%..........................................................................
% Note:
% We subdivide each sph. triangle in the list LIST 1 in 4 triangles, that
% we call "sons".
%
% Each initial planar triangle defining the generical sph. triangle, has
% 3 midpoints. In view of the midpoints, we get 4 triangles.
%
% Of each son triangle we make a structure:
%
% * L2_data.vertices : matrix 3 x 3, whose l-th row are the
%           cartesian coordinates of the l-th vertex of the j-th sub sph.
%           triangle;
% * L2_data.area     : scalar representing the area of the j-th son sph.
%           triangle;
% * L2_data.integral : scalar representing the integral of "f" on the j-th
%           son sph. triangle.
%..........................................................................

L2_data_temp=generate_sphtri_sons(vertices,f);
L2_data{end+1}=L2_data_temp;

% ... Evaluate absolute error ...
for j=1:4, L2_data_integrals(j)=L2_data_temp{j}.integral; end

% better approximation of integral on triangle (more evals)
Itemp_high=sum(L2_data_integrals);

L1_errors=[L1_errors; abs(Itemp_high-Itemp_low)];








function Ihigh=Ihigh_computation(L2_data)

%--------------------------------------------------------------------------
% Object:
% Approximation of the integral in view of approximations in LIST 2.
%--------------------------------------------------------------------------
% Input:
% L2_data: struct defining LIST 2 sph. triangles vertices, areas and
%          integral approximations.
%--------------------------------------------------------------------------
% Output:
% Ihigh:   approximation of integral of "f" over the initial spherical
%          triangle.
%--------------------------------------------------------------------------

M=length(L2_data);
Ihigh=0;
for j=1:M
    L2_data_temp=L2_data{j};
    for jj=1:4
        L2_data_temp_sons_integral(jj)=L2_data_temp{jj}.integral;
    end
    Ihigh=Ihigh+sum(L2_data_temp_sons_integral);
end











function L2_data=generate_sphtri_sons(vertices,f)

%--------------------------------------------------------------------------
% Input:
% vertices : it is a 3 x 3 matrix, where the k-th row represents the
%            cartesian coordinates of the k-th vertex (counterclockwise);
% f        : function to integrate over the spherical triangle defined by
%            vertices;
%--------------------------------------------------------------------------
% Output:
% L2_data :
%--------------------------------------------------------------------------

OA=(vertices(1,:)); OB=(vertices(2,:)); OC=(vertices(3,:));
R=norm(OA);

% ........................ compute midpoints ..............................

OAB_mid=(OA+OB)/2; OAB_mid=R*OAB_mid/norm(OAB_mid);
OAC_mid=(OA+OC)/2; OAC_mid=R*OAC_mid/norm(OAC_mid);
OBC_mid=(OB+OC)/2; OBC_mid=R*OBC_mid/norm(OBC_mid);

% ........................ triangle data ..................................

vertices=[OA; OAB_mid; OAC_mid];

L2_data{1}=make_L2_data(vertices,f);
L2_data{2}=make_L2_data([OAB_mid; OB; OBC_mid],f);
L2_data{3}=make_L2_data([OBC_mid; OC; OAC_mid],f);
L2_data{4}=make_L2_data([OAB_mid; OBC_mid; OAC_mid],f);









function L2_dataL=make_L2_data(vertices,f)

%--------------------------------------------------------------------------
% Input:
% vertices : it is a 3 x 3 matrix, where the k-th row represents the
%            cartesian coordinates of the k-th vertex (counterclockwise);
% f        : function to integrate over the spherical triangle defined by
%            vertices.
%--------------------------------------------------------------------------
% Output:
% L2_dataL : determine from the spherical triangle defined by vertices, a
%            struct data, with form:
%            * L2_dataL.vertices,
%            * L2_dataL.area,
%            * L2_dataL.integral.
%--------------------------------------------------------------------------

% vertices
L2_dataL.vertices=vertices;
P1=(vertices(1,:))'; P2=(vertices(2,:))'; P3=(vertices(3,:))';
% xyzw = sphtriquad(4,5,P1,P2,P3,0);
xyzw= cub_sphtri(4,P1,P2,P3);
L2_dataL.integral=(xyzw(:,4))'*feval(f,xyzw(:,1),xyzw(:,2),xyzw(:,3));










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


% .................... compute barycenter vertices ........................

R=norm(vertices(1,:)); vertices=vertices/R;
if  anynan(vertices(:,1))
    id = 1:size(vertices,1); idx = isnan(vertices(:,1)); nanidx = id(idx');
    CC=mean(vertices(1:(nanidx-1),:),1); CC=CC/norm(CC);
else
    CC=mean(vertices,1); CC=CC/norm(CC);
end


% .................... rotation matrix to the north pole ..................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1];
rotmat=R1*R2; inv_rotmat=rotmat';

% ........................ rotate vertices ................................

vertices_NP=(rotmat*vertices')';

% ....... map radially vertices to plane tangent to north pole ............

XX_NP=vertices_NP(:,1); YY_NP=vertices_NP(:,2); ZZ_NP=vertices_NP(:,3);
rat=1./(1+ZZ_NP);

XX_NPm=rat.*XX_NP; YY_NPm=rat.*YY_NP; ZZ_NPm=ones(size(XX_NPm));

% ............ polyshape of projected polygon to North Pole ...............

PG = polyshape(XX_NPm,YY_NPm);
PG = simplify(PG);

% ............ triangulation of projected polygon to North Pole ...........

tri = triangulation(PG);

tri_vertices_NPm=tri.Points;

phi = @(x,y) [2*x./(1+x.^2+y.^2), 2*y./(1+x.^2+y.^2), (1-x.^2-y.^2)./(1+x.^2+y.^2)];
tri_vertices_NPm=phi(tri_vertices_NPm(:,1),tri_vertices_NPm(:,2));

% ......................... output data ...................................

%tri_vertices_NP=tri_vertices_NPm./rads;
tri_vertices=R*(inv_rotmat*tri_vertices_NPm')';
tri_conn_list=tri.ConnectivityList;
