
function xyzw = cub_sphtri(n,P1,P2,P3,pos)

%--------------------------------------------------------------------------
% Object:
% The routine computes a cubature formula "xyzw" with (numerical) algebraic 
% degree of precision "n" on the spherical triangle with vertices "P1",  
% "P2", "P3", by a 2D formula of degree of precision "m+n" on the xy  
% projection of the spherical triangle rotated to put the barycenter at the
% north pole.
%
% The parameter "m" depends on the size of the circumradius and it is "2" 
% for "small" "r", becoming bigger as "r" approaches "1".
%--------------------------------------------------------------------------
% Input: 
% n: algebraic degree of precision of the rule;
% P1,P2,P3: column arrays of the spherical triangle vertices coords
% pos: positive weights for pos=1, possible neg weights (faster) for pos=0
%--------------------------------------------------------------------------
% Output:
% xywz : 4-column array of nodes cartesian coords and weights
%--------------------------------------------------------------------------
% Routines called:
% 1. compute_m (attached)
% 2. cub_circsect (external)
% 3. dCATCH (external)
%--------------------------------------------------------------------------
% Data:
% The original routine has been written by M. Vianello on May 2019.
% Last Update: 01/01/2021 by A. Sommariva.
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
%% Date: January 01, 2021
%--------------------------------------------------------------------------


% ............................ troubleshooting ............................

% working with vertices as column vectors
if size(P1,1) == 1, P1=P1'; end
if size(P2,1) == 1, P2=P2'; end
if size(P3,1) == 1, P3=P3'; end

% ------------------------------- main code -------------------------------
% We scale the problem in the unit-sphere with center [0 0 0].
% In other words, we make up a cubature rule on the unit sphere (by
% contraction w.r.t. the unit sphere and map back the so obtained rule to
% the starting sphere with radius R).

% ... barycenter ...
RR=norm(P1);
vert0=[P1 P2 P3]/RR; 
CC=1/3*(P1+P2+P3); CC=CC/norm(CC);

% .................... rotation matrix to the north pole ..................

[az,el,r] = cart2sph(CC(1),CC(2),CC(3));
phi=az; theta=pi/2-el;
cp=cos(phi); sp=sin(phi); ct=cos(theta); st=sin(theta);
R1=[ct 0 -st; 0 1 0; st 0 ct]; R2=[cp sp 0; -sp cp 0; 0 0 1]; 
R=R1*R2; invR=R';

% ............... vertices of the triangle at the North Pole ..............

vert1=R*vert0; 

% ...... computing "m" needed in determining the degree of precision ......
m=compute_m(vert1');

% .................... determining quadrature rule ........................  

xyzw=[]; % nodes on the sphere
vert1=[vert1 vert1(:,1)];

for i=1:3
    
    % affine transformation matrix 
    P=vert1(:,i); Q=vert1(:,i+1);
    om=acos(P'*Q);
    xi=cos(om/2); eta=sin(om/2);
    M=[xi eta 0 0; 0 0 xi eta; xi -eta 0 0; 0 0 xi -eta];
    h=[P(1); P(2); Q(1); Q(2)];
    u=M\h;
    T=[u(1) u(2); u(3) u(4)];
    
    % nodes and weights for the xy-projection of the rotated spher.triangle
    nw=cub_circsect(m+n,om/2,0,1); nod2=T*nw(:,1:2)';
    
    % spherical triangle nodes and weights on North Pole spher.triangle
    x=nod2(1,:); y=nod2(2,:); z=sqrt(1-x.^2-y.^2); 
    
    weights=abs(det(T))*nw(:,3)./z';
    
    % inverse rotation to original spherical triangle
    nodes=invR*[x; y; z];
    
    % cubature rule update
    xyzw=[xyzw;[nodes' weights]];
    
end

% Exporting results to the sphere with radius "R".
X=RR*xyzw(:,1); Y=RR*xyzw(:,2); Z=RR*xyzw(:,3); W=RR*xyzw(:,4);
xyzw=[X Y Z W];






function m=compute_m(vert)

%--------------------------------------------------------------------------
% Object:
% It computes the "m" positive integer value depending from the sph. tri.
% with vertices "vert", so that a WAM over the sph. triangle with degree
% equal to "n" can be obtained via WAM of degree "m+n" on its projection on
% the xy-plane.
%--------------------------------------------------------------------------
% Input:
% vert: points defining the sph. triangle (the k-th point is described 
%    by the k-th row.
%--------------------------------------------------------------------------
% Output:
% m: positive integer value depending from the sph.triangle with vertices 
% "vert", so that a WAM over the sph. triangle with degree equal to "n" can
% be obtained via cubature of degree "m+n" on its projection on the 
% xy-plane.
%--------------------------------------------------------------------------
% Data:
% First version: 13/01/2021 by A. Sommariva and M. Vianello.
%--------------------------------------------------------------------------

r=sqrt( (vert(:,1)).^2 + (vert(:,2)).^2 ); 
r0=max(r);

intvf=[0,r0]; F=@(r) sqrt(1-r); f=chebfun(F,intvf);
m=max(1,(length(f)-1));






function tw=trigauss(n,alpha,beta,method)

%--------------------------------------------------------------------------
% AUTHORS.
%--------------------------------------------------------------------------
% Alvise Sommariva and Marco Vianello, University of Padova
% July 25, 2016.
%
% previous versions with the help of G. Da Fies.
%--------------------------------------------------------------------------
% OBJECT.
%--------------------------------------------------------------------------
% It computes the angles and weights of a trigonometric quadrature formula
% on [alpha,beta], 0<beta-alpha<=pi, matching the moments up to 10^(-14).
%
% Depending on the choice of the variable 'method', some methods are
% implemented.
%
% The formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300.
%--------------------------------------------------------------------------
% INPUT.
%--------------------------------------------------------------------------
% n: trigonometric degree of exactness.
%
% [alpha,beta]: angular interval, 0<(beta-alpha)/2<=pi.
%
% method:
%         1. 'classic' implements classical trigauss rule.
%         2. 'legendre' implements (shifted) Gauss-Legendre rule
%             matching the trigonometric moments up to 10^(-14). Few nodes
%             if angles are small otherwise the cardinality might be higher
%             than in 'classic'.
%         3. 'better': chooses 'classic' or 'legendre', so to have the
%             minimal number of points, still matching the moments up to
%             10^(-14).
%         4.  'antigauss': antigaussian formula based on Legendre rule of  
%             degree M=M(n). 
%             It provides a (M+1) x 4 matrix whose entries are 
%                             tw=[tAGL wAGL tGLe wGLe]
%             where 
%             a) (tAGL,wAGL) is the antigaussian rule;
%             b) (tGLe(1:end-1),wGLe(1:end-1)) is the gaussian rule;
%         5.  'kronrod': Gauss-Kronrod formula based on Legendre rule of  
%             degree M=M(n). 
%             It provides a (2*M+1) x 3 matrix whose entries are 
%                             tw=[tGK wGK wGLf]
%             where 
%             a) (tGK,wGK) is the Gauss-Kronrod rule;
%             b) (tGK(2:2:end),wGLf(2:2:end)) is the gaussian rule;
%         6.  'antitrigauss': antigaussian formula based on trigonometric  
%             gaussian rule of degree n (whose cardinality is n+1). 
%             It provides a (n+2) x 4 matrix whose entries are 
%                             tw=[tAGL wAGL tGLe wGLe];
%             where 
%             a) (tAGL,wAGL) is the antigaussian rule;
%             b) (tGLe(1:end-1),wGLe(1:end-1)) is the gaussian rule;
%         7.  'trig_kronrod': Gauss-Kronrod based on trigonometric gaussian 
%             rule of degree n. 
%             It provides a (2*(n+1)+1) x 3 matrix whose entries are 
%                             tw=[tTK wTK wTf];
%             where 
%             a) (tTK,wTK) is the Gauss-Kronrod rule;
%             b) (tTK(2:2:end),wTf(2:2:end)) is the trig. gaussian rule;
%         otherwise: if no string is given it chooses 'better' option by 
%             default.
%--------------------------------------------------------------------------
% OUTPUT.
%--------------------------------------------------------------------------
% tw: 1) for 'classic', 'legendre', 'better', it is a 
%     M x 2 array of (angles,weights) (M depends on the choosen rule)
%     2) for 'antigauss' and 'antitrigauss' it is a M x 4 matrix (see
%     explanation above at points 4. and 6.)
%     3) for 'kronrod' and 'trig_kronrod' it is a M x 3 matrix (see
%     explanation above at points 5. and 7.)
%--------------------------------------------------------------------------
% EXAMPLES.
%--------------------------------------------------------------------------
%
% >> format long e;
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'better');
% >> twB
% 
% twB =
% 
%     -1.508419779244729e-02     1.590094129717467e-03
%     -1.251400776401954e-02     3.493153120682099e-03
%     -8.255043791082401e-03     4.927692470361337e-03
%     -2.881384626391015e-03     5.697023547188071e-03
%      2.881384626391019e-03     5.697023547188069e-03
%      8.255043791082396e-03     4.927692470361324e-03
%      1.251400776401955e-02     3.493153120682099e-03
%      1.508419779244729e-02     1.590094129717463e-03
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'classic');
% >> twB
% 
% twB =
% 
%     -1.536597391335750e-02     8.744544324142629e-04
%     -1.393392018525980e-02     1.972635825710960e-03
%     -1.146915301262857e-02     2.926255469814080e-03
%     -8.153889677130131e-03     3.662992813352760e-03
%     -4.233938891704753e-03     4.128095270647346e-03
%      2.178667008161081e-18     4.287058912019128e-03
%      4.233938891704759e-03     4.128095270647326e-03
%      8.153889677130132e-03     3.662992813352761e-03
%      1.146915301262857e-02     2.926255469814094e-03
%      1.393392018525980e-02     1.972635825710951e-03
%      1.536597391335750e-02     8.744544324142609e-04
% 
% >> format short e
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'antigauss');
% >> twB
% 
% twB =
% 
%   -1.5612e-02   5.4060e-04  -1.5084e-02   1.5901e-03
%   -1.4036e-02   2.5826e-03  -1.2514e-02   3.4932e-03
%   -1.0564e-02   4.2828e-03  -8.2550e-03   4.9277e-03
%   -5.6645e-03   5.4042e-03  -2.8814e-03   5.6970e-03
%   -7.0636e-19   5.7956e-03   2.8814e-03   5.6970e-03
%    5.6645e-03   5.4042e-03   8.2550e-03   4.9277e-03
%    1.0564e-02   4.2828e-03   1.2514e-02   3.4932e-03
%    1.4036e-02   2.5826e-03   1.5084e-02   1.5901e-03
%    1.5612e-02   5.4060e-04            0            0
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'kronrod');
% >> twB
% 
% twB =
% 
%   -1.5604e-02   2.7995e-04            0
%   -1.5084e-02   7.7659e-04   1.5901e-03
%   -1.4045e-02   1.2956e-03            0
%   -1.2514e-02   1.7537e-03   3.4932e-03
%   -1.0561e-02   2.1404e-03            0
%   -8.2550e-03   2.4607e-03   4.9277e-03
%   -5.6659e-03   2.7029e-03            0
%   -2.8814e-03   2.8494e-03   5.6970e-03
%    3.4551e-18   2.8973e-03            0
%    2.8814e-03   2.8494e-03   5.6970e-03
%    5.6659e-03   2.7029e-03            0
%    8.2550e-03   2.4607e-03   4.9277e-03
%    1.0561e-02   2.1404e-03            0
%    1.2514e-02   1.7537e-03   3.4932e-03
%    1.4045e-02   1.2956e-03            0
%    1.5084e-02   7.7659e-04   1.5901e-03
%    1.5604e-02   2.7995e-04            0
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'antitrigauss');
% >> twB
% 
% twB =
% 
%   -1.5366e-02   8.7445e-04  -1.5298e-02   1.0473e-03
%   -1.3934e-02   1.9726e-03  -1.3588e-02   2.3476e-03
%   -1.1469e-02   2.9263e-03  -1.0672e-02   3.4414e-03
%   -8.1539e-03   3.6630e-03  -6.8077e-03   4.2296e-03
%   -4.2339e-03   4.1281e-03  -2.3385e-03   4.6420e-03
%    2.1787e-18   4.2871e-03   2.3385e-03   4.6420e-03
%    4.2339e-03   4.1281e-03   6.8077e-03   4.2296e-03
%    8.1539e-03   3.6630e-03   1.0672e-02   3.4414e-03
%    1.1469e-02   2.9263e-03   1.3588e-02   2.3476e-03
%    1.3934e-02   1.9726e-03   1.5298e-02   1.0473e-03
%    1.5366e-02   8.7445e-04            0            0
% 
% >> a=-pi/200; b=pi/200; n=10; twB=trigauss(n,a,b,'trig_kronrod');
% >> twB
% 
% twB =
% 
%   -1.5640e-02   1.8370e-04            0
%   -1.5298e-02   5.1143e-04   1.0473e-03
%   -1.4611e-02   8.6012e-04            0
%   -1.3588e-02   1.1787e-03   2.3476e-03
%   -1.2265e-02   1.4628e-03            0
%   -1.0672e-02   1.7183e-03   3.4414e-03
%   -8.8397e-03   1.9398e-03            0
%   -6.8077e-03   2.1160e-03   4.2296e-03
%   -4.6243e-03   2.2427e-03            0
%   -2.3385e-03   2.3207e-03   4.6420e-03
%   -4.7437e-18   2.3475e-03            0
%    2.3385e-03   2.3207e-03   4.6420e-03
%    4.6243e-03   2.2427e-03            0
%    6.8077e-03   2.1160e-03   4.2296e-03
%    8.8397e-03   1.9398e-03            0
%    1.0672e-02   1.7183e-03   3.4414e-03
%    1.2265e-02   1.4628e-03            0
%    1.3588e-02   1.1787e-03   2.3476e-03
%    1.4611e-02   8.6012e-04            0
%    1.5298e-02   5.1143e-04   1.0473e-03
%    1.5640e-02   1.8370e-04            0
%
%--------------------------------------------------------------------------
% NOTE.
%--------------------------------------------------------------------------
% For examples about the usage of formulas of antigaussian or Kronrod type,
% see the file "demo_trigauss_error.m".
%--------------------------------------------------------------------------
%% Copyright (C) 2016
%% Gaspare Da Fies, Alvise Sommariva, Marco Vianello.
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
%% Gaspare Da Fies, Alvise Sommariva, Marco Vianello.
%%
%% Date: JULY 25, 2016
%--------------------------------------------------------------------------

if nargin <= 3
    method = 'classic';
end

validStrings = {'classic','legendre','better','antigauss',...
    'kronrod','antitrigauss','trig_kronrod'};

if ( any(strcmpi(method, validStrings)) )
    
    if ( strcmpi(method, 'classic') )
        tw=trigauss_classic(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'legendre') )
        if beta-alpha==2*pi
            tw=circle_trigquad(n,alpha,beta);
        else
            NN=degree_trigGL(n,alpha,beta);
            tw=gauss_legendre(NN,alpha,beta);
        end
    end
    
    if ( strcmpi(method, 'better') )
        if beta-alpha==2*pi
            tw=circle_trigquad(n,alpha,beta);
        else
            NN=degree_trigGL(n,alpha,beta);
            if NN <= n
                tw=gauss_legendre(NN,alpha,beta);
            else
                tw=trigauss_classic(n,alpha,beta);
            end
        end
    end
    
    if ( strcmpi(method, 'antigauss') )
        NN=degree_trigGL(n,alpha,beta);
        tw=antigauss_full(NN,alpha,beta);
    end

    if ( strcmpi(method, 'kronrod') )
        NN=degree_trigGL(n,alpha,beta);
        tw=gauss_kronrod_full(NN,alpha,beta);
    end
    
    if ( strcmpi(method, 'antitrigauss') )
        tw=antitrigauss_full(n,alpha,beta);
    end
    
    if ( strcmpi(method, 'trig_kronrod') )
        tw=trigauss_kronrod_full(n,alpha,beta);
    end
    
else
    warning('wrong string as method, choosen classical trigauss');
    tw=trigauss(n,alpha,beta,'classic');
end





function NN=degree_trigGL(n,alpha,beta)

% this routine computes the number of nodes that Gauss-Legendre needs so to
% match the trigonometric moments up to 10^(-14).

theta=(beta-alpha)/2;

if n < 500
    
    u=[1 (25:25:500)];
    
    L=[8 29    45    60    75    89   103   119   132   146   158   173  ...
        185   199 212   225   240   253   264   279   291];
    
    s=spline(u,L);
    
    NN=ceil(ppval(s,n*theta));
    
else
    
    u=n*theta;
    s=0.54*u+21;
    NN=ceil(s);
    
end



function tw=gauss_kronrod_full(n,alpha,beta)
n0=ceil(3*n/2)+1;
ab0=r_jacobi(n0,0,0);

% KRONROD.
xw=kronrod(n,ab0);
xGK=xw(:,1); wGK=xw(:,2);
tGK=(beta+alpha)/2+(beta-alpha)*xGK/2;
wGK=(beta-alpha)*wGK/2;


% GAUSS.
xwGL=gauss(n,ab0);
xGL=xwGL(:,1); wGL=xwGL(:,2);
tGL=(beta+alpha)/2+(beta-alpha)*xGL/2;

wGLl=(beta-alpha)*wGL/2; 
wGLf=zeros(size(wGK));

wGLf(2:2:end)=wGLl;

tw=[tGK wGK wGLf];





function tw=antigauss_full(n,alpha,beta)

ab=r_jacobi(n+1,0,0);
ab(end,2)=2*ab(end,2);
xwAGL=gauss(n+1,ab);

xAGL=xwAGL(:,1); wAGL=xwAGL(:,2);
tAGL=(beta+alpha)/2+(beta-alpha)*xAGL/2;
wAGL=(beta-alpha)*wAGL/2;

xwGL=gauss(n,ab(1:end-1,:));
xGL=xwGL(:,1); wGL=xwGL(:,2);
tGL=(beta+alpha)/2+(beta-alpha)*xGL/2; tGLe=[tGL; 0];
wGL=(beta-alpha)*wGL/2; wGLe=[wGL; 0];

tw=[tAGL wAGL tGLe wGLe];







function tw=trigauss_kronrod_full(n,alpha,beta)

omega=(beta-alpha)/2;
n0=ceil(3*n/2)+1;
ab0=r_trigauss(n0,omega);

% KRONROD.
xw=kronrod(n,ab0);
xTK=xw(:,1); wTK=xw(:,2);
tTK(:,1)=2*asin(sin(omega/2)*xTK(:,1))+(beta+alpha)/2;


% GAUSS.
xwT=gauss(n,ab0);
xT=xwT(:,1); wT=xwT(:,2);
tT(:,1)=2*asin(sin(omega/2)*xT(:,1))+(beta+alpha)/2;

wTf=zeros(size(wTK));

wTf(2:2:end)=wT;

tw=[tTK wTK wTf];





function tw=antitrigauss_full(n,alpha,beta)

omega=(beta-alpha)/2;

ab=r_trigauss(n+1,omega);
ab(end,2)=2*ab(end,2);
xwAT=gauss(n+1,ab);

xAT=xwAT(:,1); wAT=xwAT(:,2);
tAT=2*asin(sin(omega/2)*xwAT(:,1))+(beta+alpha)/2;

xwT=gauss(n,ab(1:end-1,:));
xT=xwT(:,1); wT=xwT(:,2);
tT=2*asin(sin(omega/2)*xT)+(beta+alpha)/2;
tTe=[tT; 0];
wTe=[wT; 0];

tw=[tAT wAT tTe wTe];




function tw=gauss_legendre(n,alpha,beta)

% [xGL, wGL] = legpts(N); wGL=wGL';
ab=r_jacobi(n,0,0);
xw=gauss(n,ab); xGL=xw(:,1); wGL=xw(:,2);
t=(beta+alpha)/2+(beta-alpha)*xGL/2;
w=(beta-alpha)*wGL/2;
tw=[t w];




function tw=trigauss_classic(n,alpha,beta)

% by Gaspare Da Fies and Marco Vianello, University of Padova
% 8 Nov 2011

% computes the n+1 angles and weights of a trigonometric gaussian
% quadrature formula on [alpha,beta], 0<beta-alpha<=pi

% uses the routines chebyshev.m, gauss.m from
% www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html
% we suggest to put the following statements
% ab = zeros(N,2); sig = zeros(N+1,2*N);
% at the beginning of the body of chebyshev.m to speed-up execution

% input:
% n: trigonometric degree of exactness
% [alpha,beta]: angular interval, 0<beta-alpha<=pi

% output:
% tw: (n+1) x 2 array of (angles,weights)

% the formula integrates the canonical trigonometric basis with accuracy
% from about 10^(-15) (small omega) to about 10^(-13) (omega-->pi)
% up to n=300


% half-length of the angular interval
omega=(beta-alpha)/2;

if omega == pi
    tw=circle_trigquad(n,0,2*pi);
else
    
    ab=r_trigauss(n,omega);
    
    % Gaussian formula for the weight function above
    xw=gauss(n+1,ab);
    
    % angles and weights for the trigonometric gaussian formula
    tw(:,1)=2*asin(sin(omega/2)*xw(:,1))+(beta+alpha)/2;
    tw(:,2)=xw(:,2);
    
end






function ab=r_trigauss(n,omega)

% modified Chebyshev moments by recurrence
z(1)=2*omega;

z(n+1)=quadgk(@(t)cos(2*n*acos(sin(t/2)/sin(omega/2))),...
    -omega,omega,'MaxIntervalCount',50000);
temp=(2:2:2*n-1);
dl=1/4-1./(4*(temp-1));
dc=1/2-1/sin(omega/2)^2-1./(2*(temp.^2-1));
du=1/4+1./(4*(temp+1));
d=4*cos(omega/2)/sin(omega/2)./(temp.^2-1)';
d(n-1)=d(n-1)-du(n-1)*z(n+1);
z(2:n)=tridisolve(dl(2:n-1),dc(1:n-1),du(1:n-2),d(1:n-1));
mom=zeros(1,2*n+2);
mom(1:2:2*n+1)=z(1:n+1);

% normalization of the moments (monic polynomials)
k=(3:length(mom));
mom(3:end)=exp((2-k)*log(2)).*mom(3:end);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*n+1,1);
abm(:,2)=0.25*ones(2*n+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
% w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm
[ab,normsq]=chebyshev(n+1,mom,abm);






function x = tridisolve(a,b,c,d)
%   TRIDISOLVE  Solve tridiagonal system of equations.
% From Cleve Moler's Matlab suite
% http://www.mathworks.it/moler/ncmfilelist.html

%     x = TRIDISOLVE(a,b,c,d) solves the system of linear equations
%     b(1)*x(1) + c(1)*x(2) = d(1),
%     a(j-1)*x(j-1) + b(j)*x(j) + c(j)*x(j+1) = d(j), j = 2:n-1,
%     a(n-1)*x(n-1) + b(n)*x(n) = d(n).
%
%   The algorithm does not use pivoting, so the results might
%   be inaccurate if abs(b) is much smaller than abs(a)+abs(c).
%   More robust, but slower, alternatives with pivoting are:
%     x = T\d where T = diag(a,-1) + diag(b,0) + diag(c,1)
%     x = S\d where S = spdiags([[a; 0] b [0; c]],[-1 0 1],n,n)

x = d;
n = length(x);
for j = 1:n-1
    mu = a(j)/b(j);
    b(j+1) = b(j+1) - mu*c(j);
    x(j+1) = x(j+1) - mu*x(j);
end
x(n) = x(n)/b(n);
for j = n-1:-1:1
    x(j) = (x(j)-c(j)*x(j+1))/b(j);
end












% CHEBYSHEV Modified Chebyshev algorithm.
%
%    Given a weight function w encoded by its first 2n modified
%    moments, stored in the (row) vector mom, relative to monic 
%    polynomials defined by the (2n-1)x2 array abm of their
%    recurrence coefficients, [ab,normsq]=CHEBYSHEV(n,mom,abm)
%    generates the array ab of the first n recurrence coefficients
%    of the orthogonal polynomials for the weight function w, and 
%    the vector normsq of their squared norms. The n alpha-
%    coefficients are stored in the first column, the n beta-
%    coefficients in the second column, of the nx2 array ab. The
%    call [ab,normsq]=CHEBYSHEV(n,mom) does the same, but using the 
%    classical Chebyshev algorithm. If n is larger than the sizes
%    of mom and abm warrant, then n is reduced accordingly.
%
function [ab,normsq]=chebyshev(N,mom,abm)
if N<=0, error('N out of range'), end
if N>size(mom,2)/2, N=size(mom,2)/2; end
if nargin<3, abm=zeros(2*N-1,2); end
if N>(size(abm,1)+1)/2, N=(size(abm,1)+1)/2; end
ab(1,1)=abm(1,1)+mom(2)/mom(1); ab(1,2)=mom(1);
if N==1, normsq(1)=mom(1); return, end
sig(1,1:2*N)=0; sig(2,:)=mom(1:2*N);
for n=3:N+1
  for m=n-1:2*N-n+2
    sig(n,m)=sig(n-1,m+1)-(ab(n-2,1)-abm(m,1))*sig(n-1,m) ...
      -ab(n-2,2)*sig(n-2,m)+abm(m,2)*sig(n-1,m-1);
  end
  ab(n-1,1)=abm(n-1,1)+sig(n,n)/sig(n,n-1)-sig(n-1,n-1)/ ...
    sig(n-1,n-2);
  ab(n-1,2)=sig(n,n-1)/sig(n-1,n-2);
end
for n=1:N, normsq(n)=sig(n+1,n); end; normsq=normsq'; 






function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
  J(n,n-1)=sqrt(ab(n,2));
  J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];








function xyw = cub_circsect(n,omega,r1,r2)

%--------------------------------------------------------------------------
% Object:
% The routine computes the nodes and weights of a product gaussian
% formula on a circular annular sector centered at the origin
% with angles in [-omega,omega]
%--------------------------------------------------------------------------
% Input:
% n: algebraic degree of exactness
% omega: half-length of the angular interval, 0<omega<=pi
% r1,r2: internal and external radius, 0<=r1<r2
%--------------------------------------------------------------------------
% Output:
% xyw: (ceil((n+2)/2) x (n+1)) x 3 array of (xnodes,ynodes,weights)
%--------------------------------------------------------------------------
% Required routines:
% 1. r_jacobi.m (www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html)
% 2. gauss.m (www.cs.purdue.edu/archives/2002/wxg/codes/OPQ.html)
% 3. quad_trig.m 
%--------------------------------------------------------------------------
% Written by Gaspare Da Fies and Marco Vianello, University of Padova
% Date: November 8, 2011.
% Last update: January 4. 2020.
%--------------------------------------------------------------------------

% trigonometric gaussian formula on the arc
tw=quad_trig(n,-omega,omega);

% algebraic gaussian formula on the radial segments
ab=r_jacobi(ceil((n+2)/2),0,0);
xw=gauss(ceil((n+2)/2),ab);
xw(:,1)=xw(:,1)*(r2-r1)/2+(r2+r1)/2;
xw(:,2)=xw(:,2)*(r2-r1)/2;

% creating the polar grid
[r,theta]=meshgrid(xw(:,1),tw(:,1));
[w1,w2]=meshgrid(xw(:,2),tw(:,2));

% nodal cartesian coordinates and weights
xyw(:,1)=r(:).*cos(theta(:));
xyw(:,2)=r(:).*sin(theta(:));
xyw(:,3)=r(:).*w1(:).*w2(:);







function ab=r_jacobi(N,a,b)

nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1
 ab=[nu mu]; return
end

N=N-1;
n=1:N;
nab=2*n+a+b;
nuadd=(b^2-a^2)*ones(1,N)./(nab.*(nab+2));
A=[nu nuadd];
n=2:N;
nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
abadd=[mu; B1; B'];
ab=[A' abadd];






function tw=quad_trig(n,alpha,beta)

%--------------------------------------------------------------------------
% Object:
% Computation of trigonometric gaussian rules on the unit arc from "alpha"
% to "beta".
%--------------------------------------------------------------------------
% Inputs:
% n    : the rule computes computes the n+1 angles and weights of a 
%        trigonometric gaussian quadrature formula on [alpha,beta], with 
%                          0 < beta-alpha <= 2*pi;
% alpha, beta: arc angles for trigonometric gaussian rules on the unit arc
%       from "alpha" to "beta".
%--------------------------------------------------------------------------
% Outputs
% tw   : (n+1) x 2 matrix, where the first column contains the nodes, while
%        the second one contains the weights.
%--------------------------------------------------------------------------
% First version: May 18, 2013 by G. Da Fies, A. Sommariva, M. Vianello.
% Successive versions have been made by G. Meurant, A. Sommariva and
% M. Vianello.
% Last release: January 04, 2020.
%--------------------------------------------------------------------------

% .......................... troubleshooting  .............................

if nargin < 1, n=10; end
if nargin < 2, beta=pi; alpha=-beta; end
if nargin < 3, beta=alpha; alpha=-beta; end

n=n+1;
omega=(beta-alpha)/2;
ab = r_subchebyshev(n,omega);
xw_symm_eigw = SymmMw(n,ab);
tw=quad_trig_conversion(xw_symm_eigw,omega);
tw(:,1)=tw(:,1)+(beta+alpha)/2;









function ab=r_subchebyshev(n,omega)

%--------------------------------------------------------------------------
% Object:
% Recurrence coeffs for the monic OPS w.r.t. the weight function
%         w(x)=2*sin(omega/2)/sqrt(1-sin^2(omega/2)*x^2)
% by the modified Chebyshev algorithm.
% The reference angle of the rule is [-omega,omega].
%--------------------------------------------------------------------------
% Inputs:
% n     : number of points.
% omega : arc angle.
%--------------------------------------------------------------------------
% Output:
% ab   : three terms recursion
%--------------------------------------------------------------------------
% First version: May 18, 2013 by G. Da Fies, A. Sommariva, M. Vianello.
% Successive versions have been made by G. Meurant, A. Sommariva and
% M. Vianello.
% Last release: January 04, 2020.
%--------------------------------------------------------------------------

N=n; n=n-1;

% modified Chebyshev moments by recurrence
if rem(N,2) == 1, NN=N+1; nn=n+1; else NN=N; nn=n; end
mom=fast_moments_computation(omega,2*nn+1);

% recurrence coeffs of the monic Chebyshev polynomials
abm(:,1)=zeros(2*nn+1,1);
abm(:,2)=0.25*ones(2*nn+1,1); abm(1,2)=pi; abm(2,2)=0.5;

% recurrence coeffs for the monic OPS w.r.t. the weight function
ab = fast_chebyshev(NN,mom,abm);











function ab=fast_chebyshev(N,mom,abm)

%--------------------------------------------------------------------------
% Object:
% Modified Chebyshev algorithm, that works only for the subperiodic weight 
% function.
%--------------------------------------------------------------------------
% From Gautschi's code (simplified)
% Mar 2012
%--------------------------------------------------------------------------
% Optimized version by G. Meurant.
%--------------------------------------------------------------------------

ab = zeros(N,2);
sig = zeros(N+1,2*N);

ab(1,2) = mom(1);

sig(1,1:2*N) = 0;
sig(2,:) = mom(1:2*N);

for n = 3:N+1
    for m = n-1:2*N-n+2
        sig(n,m) = sig(n-1,m+1) + abm(m,2) * sig(n-1,m-1) - ...
            ab(n-2,2) * sig(n-2,m);
    end
    
    ab(n-1,2) = sig(n,n-1) / sig(n-1,n-2);
end









function mom=fast_moments_computation(omega,n)

%--------------------------------------------------------------------------
% Object: 
%--------------------------------------------------------------------------
% Inputs:
%--------------------------------------------------------------------------
% Outputs:
%--------------------------------------------------------------------------
% Authors G. Meurant and A. Sommariva
% June 2012
%--------------------------------------------------------------------------

mom=zeros(1,n+1);
mom(1)=2*omega; % FIRST MOMENT.

if(n>=2)
    
    if(omega<=1/4*pi)
        l=10;
    elseif(omega<=1/2*pi)
        l=20;
    elseif(omega<=3/4*pi)
        l=40;
    else
        if omega == pi
            l=2*ceil(10*pi);
        else
            l=2*ceil(10*pi/(pi-omega));
        end
    end
    
    
    temp=(2:2:n+2*l-2); % AUXILIAR VECTORS.
    temp2=temp.^2-1;
    
    dl=1/4 -1./(4*(temp-1)); % DIAGONALS.
    dc=1/2 -1/sin(omega/2)^2 -1./(2*temp2);
    du=1/4 +1./(4*(temp+1));
    
    d=4*cos(omega/2)/sin(omega/2)./temp2'; % COMPUTING KNOWN TERM.
    d(end)=d(end);                         % PUT LAST MOMENT NULL
    
    z=tridisolve(dl(2:end),dc,du(1:end-1),d); % SOLVE SYSTEM.
    mom(3:2:n+1)=z(1:floor(n/2)); % SET ODD MOMENTS.
    
end

mom=mom';

normalized = 0;

if normalized == 0
    M=length(mom);
    kk=2.^(-((1:2:M)-2))'; kk(1)=1;
    v=ones(M,1);
    v(1:2:M)=kk;
    mom=v.*mom;
end









function xw=SymmMw(N,ab)

%--------------------------------------------------------------------------
% Object:
% Computation of the nodes and weights for a symmetric weight function
% this version uses the reduced matrix and eig and computation of weights 
% with the 3-term recurrence.
%--------------------------------------------------------------------------
% Input:
% N : cardinality of the rule
% ab: 3-term recurrence for the orthogonal polynomials same as in OPQ
%     ab(1,2) is the 0th moment.
%--------------------------------------------------------------------------
% Output:
% xw : xw(:,1) nodes, xw(:,2) weights of the quadrature rule
%--------------------------------------------------------------------------
% Reference paper:
% Fast variants of the Golub and Welsch algorithm for symmetric
% weight functions by G. Meurant and A. Sommariva (2012)
%--------------------------------------------------------------------------
% Data:
% Written by G. Meurant and A. Sommariva on June 2012
%--------------------------------------------------------------------------

N0 = size(ab,1);
if N0 < N
    error('SymmMw: input array ab is too short')
end

na = norm(ab(:,1));
if na > 0
    error('SymmMw: the weight function must be symmetric')
end

% computation of the reduced matrix in vectors (a,b)

if mod(N,2) == 0
    even = 1;
    Nc = N / 2;
else
    even = 0;
    Nc = fix(N / 2) +1;
end


absd = ab(:,2);
absq = sqrt(absd);

a = zeros(1,Nc);
b = a;

switch even
    case 1
        % N even
        a(1) = absd(2);
        b(1) = absq(2) * absq(3);
        
        k = [2:Nc-1];
        a(k) = absd(2*k-1) + absd(2*k);
        b(k) = absq(2*k) .* absq(2*k+1);
        a(Nc) = absd(N) + absd(N-1);
        start = 1;
        
        J = diag(a) + diag(b(1:Nc-1),1) + diag(b(1:Nc-1),-1);
        t = sort(eig(J));
        w = weights_3t(t',a,b);
        % w are the squares of the first components
        w = w' / 2;
    case 0
        % N odd
        a(1) = absd(2);
        b(1) = absq(2) * absq(3);
        
        k = [2:Nc-1];
        a(k) = absd(2*k-1) + absd(2*k);
        b(k) = absq(2*k) .* absq(2*k+1);
        a(Nc) = absd(N);
        start = 2;
        
        % the first node must be zero
        J = diag(a) + diag(b(1:Nc-1),1) + diag(b(1:Nc-1),-1);
        t = sort(eig(J));
        t(1) = 0;
        w = weights_3t(t',a,b);
        w = [w(1); w(2:end)' / 2];
    otherwise
        error('this is not possible')
end

xwp = sqrt(t);

xw(:,1) = [-xwp(end:-1:start,1); xwp];
xw(:,2) = ab(1,2) * ([w(end:-1:start); w]);









function tw=quad_trig_conversion(xw,omega)

%--------------------------------------------------------------------------
% Object: 
%--------------------------------------------------------------------------
% Inputs:
%--------------------------------------------------------------------------
% Outputs:
%--------------------------------------------------------------------------
% Authors G. Meurant and A. Sommariva
% June 2012
%--------------------------------------------------------------------------

tw(:,1)=2*asin(sin(omega/2)*xw(:,1));
tw(:,2)=xw(:,2);









function w=weights_3t(t,a,b)

%--------------------------------------------------------------------------
% Object: 
% Squares of the 1st components of eigenvectors from the 3-term
% recurrence relation of the orthogonal polynomials
%--------------------------------------------------------------------------
% Inputs:
% t: nodes
% a,b: coefficients of the 3-term recurrence
%--------------------------------------------------------------------------
% Outputs
% w: squares of the first components of the eigenvectors
%--------------------------------------------------------------------------
% Authors G. Meurant and A. Sommariva
% June 2012
%--------------------------------------------------------------------------

N = length(t);

P = zeros(N,N);
P(1,:) = ones(1,N);
P(2,:) = (t - a(1)) / b(1);

for k = 3:N
    k1 = k - 1;
    k2 = k - 2;
    P(k,:) = ((t - a(k1)) .* P(k1,:) - b(k2) * P(k2,:)) / b(k1);
end

P2 = P .* P;
w = 1 ./ sum(P2);



