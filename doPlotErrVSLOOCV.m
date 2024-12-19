function doPlotErrVSLOOCV(domain_type,rbf_type, f_type) 

if nargin < 1
    domain_type = 1;
end

if nargin < 2
    rbf_type = 1;
end

if nargin< 3
    f_type = 27;
end

Nset = 1:2;
% Defining the test function
[f,~] = test_functions(f_type);

switch domain_type
    case 1
        domain = coastline_africa(0);
    case 2
        domain = coastline_australia(0);
    case 3
        domain = coastline_america;
end

Vdeg = domain.Vertices;
[Vx,Vy,Vz] = sph2cart(deg2rad(Vdeg(:,1)),deg2rad(Vdeg(:,2)),1);
vertices = [Vx,Vy,Vz];

XWC=cub_sphpgon(1,vertices);
area = sum(XWC(:,4));

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

Itrue = adaptive_cub_sphpgon(vertices,f,1e-12,1e-12);
atol=10^(-12);
rtol=atol;
avgEF = zeros(length(Nset),length(n));
avgEFCOND = NaN*ones(length(Nset),length(n));
err_rel = zeros(length(Nset),length(n));
for i=Nset
    N = i*100;
    centers = PtsSphPol(N,'H',vertices);
    DM_data = centers*centers';
    % Evaluating the centers
    fdata = f(centers(:,1),centers(:,2),centers(:,3));

    for j = 1:length(n)
        IM=phi(n(j),DM_data);
        invIM = pinv(IM);
        coeff=IM\fdata;

        EF = (coeff)./diag(invIM);

        % Compute the average of the norm of EF
        if condest(IM)<1e+18
            avgEFCOND(i,j) = abs(sum(area/size(centers,1)*EF(:)));
        end
        avgEF(i,j) = abs(sum(area/size(centers,1)*EF(:)));

        % Computing the integration of the RBF-interpolant
        I = adaptive_RBF_sphpgon(vertices,coeff,phi,n(j),centers,atol,rtol);

        err_abs = abs(I-Itrue);
        err_rel(i,j) = err_abs/abs(Itrue);
    end
end

[minavgEFCOND,idxMinavgEFCOND] = min(avgEFCOND,[],2);
[minavgEF,idxMinavgEF] = min(avgEF,[],2);

IndicesMat = sub2ind(size(avgEF),1:length(idxMinavgEF),idxMinavgEF');
IndicesMatCOND = sub2ind(size(avgEF),1:length(idxMinavgEFCOND),idxMinavgEFCOND');

figure
semilogy(100*Nset,err_rel(:,1),'b')
hold on
semilogy(Nset*100,err_rel(:,2),'r')
semilogy(Nset*100,err_rel(:,3),'m')
semilogy(Nset*100,err_rel(:,4),'k')
semilogy(Nset*100,err_rel(:,5),'Color',[0 0.4470 0.7410])
semilogy(Nset*100,err_rel(:,6),'Color',[0.8500 0.3250 0.0980])
semilogy(Nset*100,err_rel(:,7),'Color',[0.9290 0.6940 0.1250])
semilogy(Nset*100,err_rel(:,8),'Color',[0.6350 0.0780 0.1840])
semilogy(Nset*100,avgEF(:,1),'b--')
semilogy(Nset*100,avgEF(:,2),'r--')
semilogy(Nset*100,avgEF(:,3),'m--')
semilogy(Nset*100,avgEF(:,4),'k--')
semilogy(Nset*100,avgEF(:,5),'--','Color',[0 0.4470 0.7410])
semilogy(Nset*100,avgEF(:,6),'--','Color',[0.8500 0.3250 0.0980])
semilogy(Nset*100,avgEF(:,7),'--','Color',[0.9290 0.6940 0.1250])
semilogy(Nset*100,avgEF(:,8),'--','Color',[0.6350 0.0780 0.1840])
semilogy(Nset*100,minavgEF,'ko','MarkerSize',7,'MarkerFaceColor','k')
semilogy(Nset*100,err_rel(IndicesMat),'ko','MarkerSize',7,'MarkerFaceColor','k')
semilogy(Nset*100,minavgEFCOND,'ro','MarkerSize',7,'MarkerEdgeColor','r')
semilogy(Nset*100,err_rel(IndicesMatCOND),'ro','MarkerSize',7,'MarkerEdgeColor','r')
legend('1','2','3','4','5','6','7','8','Location','NorthEast')%,'FontSize',7)
axis([Nset(1)*100 Nset(end)*100 1e-10 1])
