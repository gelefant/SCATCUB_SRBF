function doPlotSphere(domain_type,f_type)

if nargin<2 f_type = 38; end
if nargin<1 domain_type = 4; end

Neval = 150000;
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

p = haltonset(2);
Xsq = net(p,N);
Xsq = [2*Xsq(:,1)-1,2*pi*Xsq(:,2)];
X = [sqrt(1-Xsq(:,1).^2).*cos(Xsq(:,2)), sqrt(1-Xsq(:,1).^2).*sin(Xsq(:,2)), Xsq(:,1)];

F = f(xeval(:,1),xeval(:,2),xeval(:,3));

bary = [mean(Vx),mean(Vy),mean(Vz)];
bary = bary/norm(bary);

figure(1)
scatter3(X(:,1),X(:,2),X(:,3),5,[244,244,244]/256,'.', 'MarkerEdgeAlpha', 0.3, 'MarkerFaceAlpha', 0.3)
hold on;
scatter3(xeval(:,1),xeval(:,2),xeval(:,3),15,F,'.')
colormap winter
% colormap(gca,map)
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
