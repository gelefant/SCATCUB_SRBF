function demo_scatcubrbf_SphAnn(testfun, domain, verbose)

if nargin<1
    testfun = 36;
end

if nargin<2
    % domain = coastline_africa(0);
    domain = holepoly;
end

if nargin<3
    verbose = 1;
end

if verbose
    fprintf('\n \n \t *** STATISTICS *** \n');
end

% Parameters
N = 100; % # scattered data in the spherical polygon
rbf_type = 2; 

% Defining the test function
[f,~] = test_functions(testfun);

% Defining the vertices of the spherical polygon
Vdeg = domain.Vertices;
V = [2*Vdeg(:,1)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2),2*Vdeg(:,2)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2),(1-Vdeg(:,1).^2-Vdeg(:,2).^2)./(1+Vdeg(:,1).^2+Vdeg(:,2).^2)];
Vx = V(:,1); Vy = V(:,2); Vz=V(:,3)+0.35;
V = [Vx,Vy,Vz]./vecnorm([Vx,Vy,Vz],2,2);
Vx = V(:,1); Vy = V(:,2); Vz = V(:,3);
vertices = [Vx,Vy,Vz];

% Defining scattered centers in the spherical cap containing the spherical
% polygon
centers = PtsSphPol(N,'H',vertices);

% Evaluating the centers
fdata = f(centers(:,1),centers(:,2),centers(:,3));

I = SCATCUB_SRBF(centers, fdata, domain, rbf_type, 1, verbose);

atol=10^(-12);
rtol=atol;

Itrue = adaptive_cub_sphpgon(vertices,f,atol,rtol);

err_abs = abs(I-Itrue);
err_rel = err_abs/abs(Itrue);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    switch rbf_type
        case 1
            rbf_string = 'RP';
        case 2
            rbf_string = 'TPS';
    end

    fprintf('\n \t RBF TYPE                   : %s',rbf_string);
    fprintf('\n \t CENTERS CARDIN.            : %4.0f',N);

    fprintf('\n \n \t FUNCTION TYPE              : %s',func2str(f));
    fprintf('\n \n  \t APPROX VALUE               : %1.15e',I);
    fprintf('\n \t EXACT VALUE                : %1.15e',Itrue);
    fprintf('\n \n \t * ABSOLUTE ERR. SCATCUB    : %1.3e',err_abs);
    if abs(I) > 0
        fprintf('\n \t * RELATIVE ERR. SCATCUB    : %1.3e',err_rel);
    end
    fprintf('\n\n')

end
