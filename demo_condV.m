function demo_condV(rbf_type, domain, do_plot, isdeg)

if nargin < 1
    rbf_type = 2;
end

if nargin < 2
    % domain = coastline_africa(0);
%     domain = coastline_australia(0);
domain_type =4;
domain = holepoly;
end

if nargin < 3
    do_plot = 1;
end

if nargin < 4
    isdeg = 0; % 1 if coordinates are in degree, 0 if are in radiant
end

% Parameters
Nset = 100:100:1600;

% Extracting the vertices of the domain
Vdeg = domain.Vertices;

if isdeg 
    Vdeg = deg2rad(Vdeg(:,1),Vdeg(:,2));
end

V = [2*Vdeg(:,1)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2),2*Vdeg(:,2)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2),(1-Vdeg(:,1).^2-Vdeg(:,2).^2)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2)];
Vx = V(:,1); Vy = V(:,2); Vz=V(:,3)+0.35;
V = [Vx,Vy,Vz]./vecnorm([Vx,Vy,Vz],2,2);
Vx = V(:,1); Vy = V(:,2); Vz = V(:,3);
vertices = [Vx,Vy,Vz];

switch rbf_type
    case 1
        phi = @(n,r) (2-2*r).^(n/2);   % Radial power (RP)
        % Parameters for RBF (power) parameter loop below
        minn = 1; maxn = 15;
        n = minn:2:maxn;
    case 2
        phi = @(n,r) (1-r).^(n).*log(2-2*r+10^(-50)); % Thin-Plate Spline (TPS)
        % Parameters for RBF (power) parameter loop below
        minn = 1; maxn = 8;
        n = minn:maxn;
end

j = 1; condV = zeros(length(n),length(Nset));
for N = Nset
% Defining scattered centers in the spherical cap containing the spherical
% polygon
centers = PtsSphPol(N,'H',vertices);

DM_data = centers*centers';

for i = 1:length(n)
    IM=phi(n(i),DM_data);
    condV(i,j) = cond(IM);
end

j = j + 1;
end

if do_plot
    figure(1)
semilogy(Nset,condV(1,:),'b')
hold on
semilogy(Nset,condV(2,:),'r')
semilogy(Nset,condV(3,:),'m')
semilogy(Nset,condV(4,:),'k')
semilogy(Nset,condV(5,:),'Color',[0 0.4470 0.7410])
semilogy(Nset,condV(6,:),'Color',[0.8500 0.3250 0.0980])
semilogy(Nset,condV(7,:),'Color',[0.9290 0.6940 0.1250])
semilogy(Nset,condV(8,:),'Color',[0.6350 0.0780 0.1840])
legend('1','2','3','4','5','6','7','8','Location','NorthWest')
axis([Nset(1) Nset(end) 1 1e30])
end