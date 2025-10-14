function doPlotErrSphere(domain_type,f_type,rbf_type)

if nargin<3 rbf_type = 2; end
if nargin<2 f_type = 38; end
if nargin<1 domain_type = 4; end

Neval = 250000;
N = 100000;

switch domain_type
    case 1
        domain = coastline_africa(0);
    case 2
        domain = coastline_australia(0);
    case 3
        domain = coastline_america;
    case 4
        domain = holepoly;
end

Vdeg = domain.Vertices;
if domain_type == 4
    V = [2*Vdeg(:,1)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2),2*Vdeg(:,2)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2),(1-Vdeg(:,1).^2-Vdeg(:,2).^2)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2)];
Vx = V(:,1); Vy = V(:,2); Vz=V(:,3)+0.35;
V = [Vx,Vy,Vz]./vecnorm([Vx,Vy,Vz],2,2);
Vx = V(:,1); Vy = V(:,2); Vz = V(:,3);
else
[Vx,Vy,Vz] = sph2cart(deg2rad(Vdeg(:,1)),deg2rad(Vdeg(:,2)),1);
end
vertices = [Vx,Vy,Vz];

[f,~]=test_functions(f_type);

xeval = PtsInSphPol(Neval,'S',vertices);

switch rbf_type
    case 1
        phi = @(n,r) (2-2*r).^(n/2);   % Radial power (RP)
        % Parameters for RBF (power) parameter loop below
        switch domain_type 
            case 1
                switch f_type
                    case 27
                        n = 3;
                    case 28
                        n = 5;
                    case 29
                        n = 5;
                    otherwise
                        n = 5;
                end
            case 2
                switch f_type
                    case 27
                        n = 5;
                    case 30
                        n = 5;
                    case 31
                        n = 5;
                    otherwise
                        n = 5;
                end
            case 3
                n = 5;
            case 4
                n = 1;
        end
    case 2
        phi = @(n,r) (1-r).^(n).*log(2-2*r+10^(-50)); % Thin-Plate Spline (TPS)
        % Parameters for RBF (power) parameter loop below
        switch domain_type 
            case 1
                switch f_type
                    case 27
                        n = 3;
                    case 28
                        n = 3;
                    case 29
                        n = 3;
                    otherwise
                        n = 3;
                end
            case 2
                switch f_type
                    case 27
                        n = 3;
                    case 30
                        n = 3;
                    case 31
                        n = 3;
                    otherwise
                        n = 3;
                end
            case 3
                n = 5;
            case 4
                n = 1;
        end
end

centers = PtsSphPol(1600,'H',vertices);
DM_data = centers*centers';
V=phi(n,DM_data);
fdata = f(centers(:,1),centers(:,2),centers(:,3));
coeff=V\fdata;

p = haltonset(2);
Xsq = net(p,N);
Xsq = [2*Xsq(:,1)-1,2*pi*Xsq(:,2)];
X = [sqrt(1-Xsq(:,1).^2).*cos(Xsq(:,2)), sqrt(1-Xsq(:,1).^2).*sin(Xsq(:,2)), Xsq(:,1)];

IM=phi(n,xeval*centers');
F = f(xeval(:,1),xeval(:,2),xeval(:,3));
I = IM*coeff;

bary = [mean(Vx),mean(Vy),mean(Vz)];
bary = bary/norm(bary);

S = F - I;

figure(1)
scatter3(X(:,1),X(:,2),X(:,3),5,[244,244,244]/256,'.', 'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)
hold on;
scatter3(xeval(:,1),xeval(:,2),xeval(:,3),15,S,'.')
colormap winter
%colormap(gca,map)
colorbar
axis equal
% plot3(bary(1),bary(2),bary(3),'.','MarkerSize', 30)
% plot3(vertices(:,1),vertices(:,2),vertices(:,3),'.','MarkerSize',12);
hold off
switch domain_type
    case 1
        view(110,15)
    case 2
        view(-130,-10)
    case 3
        view(10,10)
    case 4
        view(50,25)
end

% [I,Ierr,flag,Ihigh,iters,tri_vertices,tri_conn_list,L1_vertices]=...
%     adaptive_RBF_sphpgon(vertices,coeff,phi,n,centers,1e-12,1e-12);
% 
% figure(2)
% scatter3(xeval(:,1),xeval(:,2),xeval(:,3),15,F,'.')
% colormap winter
% hold on
% view(90, 0);
% for k=1:length(L1_vertices)
%       tri_L=L1_vertices{k};
%       tri_L(end+1,:)=tri_L(1,:);
%       plot3(tri_L(:,1),tri_L(:,2),tri_L(:,3),'r-','LineWidth',1.5);
%       hold on;
% end
