function demo_adaptive(f_type,domain_type)

if nargin<2 domain_type = 2; end

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

rtol = 1e-14;
atol = rtol;

if nargin<1
    tests = 100;
    adeV=3:3:12;

    REV=[]; logerrV=[];

    for jj=1:length(adeV)

        ade=adeV(jj);
        % Computing an cubature rule with ADE n
        XWC=cub_sphpgon(ade,vertices);

        for k=1:tests
            a=rand(1); b=rand(1); c=rand(1); d=rand(1);
            f=@(x,y,z) (a+b*x+c*y+d*z).^ade;

            [Iadapt,Iadapt_err,adapt_flag,adapt_Ihigh,adapt_iters,adapt_tri_vertices,adapt_tri_conn_list,adapt_L1_vertices]=...
                adaptive_cub_sphpgon(vertices,f,atol,rtol);

            % 3a. Integral via the rule
            fnodesC=feval(f,XWC(:,1),XWC(:,2),XWC(:,3));
            wC=XWC(:,4);
            I=wC'*fnodesC;

            RE(k,1) = abs(Iadapt-I)/abs(I);
        end

        fprintf('\n \n \t ADE   : %2.0f',ade);
        fprintf('\n \t RE max: %1.3e',max(RE));
        kpos=find(RE > 0);
        logerr=10^(sum(log10(RE(kpos)))/length(kpos));
        fprintf('\n \t RE log: %2.3e',logerr);

        REV=[REV RE];
        logerrV=[logerrV; logerr];
    end

    h=figure(1);
    f1=ishandle(h)&&strcmp(get(h,'type'),'figure'); if f1,clf(1);end
    figure(1)
    axis equal;
    for jj=1:length(adeV)
        n=adeV(jj);
        reL=REV(:,jj); log_reL=logerrV(jj);
        plot_errors(jj,n,reL,log_reL);
        hold on;
    end
    ax=gca;
    ax.XAxis.FontSize = 16;
    ax.YAxis.FontSize = 16;
    xlim([adeV(1)-1,adeV(end)+1]);
    xticks(adeV)
    hold off
    fprintf('\n \n');
else
    ade = 50;

    [f,~]=test_functions(f_type);
    test = 20;
    XWC=cub_sphpgon(ade,vertices);
    cpu_time = 0;
    for k = 1:test
    tAdapt=tic;
    [Iadapt,Iadapt_err,adapt_flag,adapt_Ihigh,adapt_iters,adapt_tri_vertices,adapt_tri_conn_list,adapt_L1_vertices]=...
        adaptive_cub_sphpgon(vertices,f,atol,rtol);
    cpu_time = cpu_time + toc(tAdapt);
    end
    cpu_time = cpu_time/test;
    % 3a. Integral via the rule
    fnodesC=feval(f,XWC(:,1),XWC(:,2),XWC(:,3));
    wC=XWC(:,4);
    I=wC'*fnodesC;
    RE = abs(Iadapt-I)/abs(I);

    fprintf('\n \t ------------------------------');
    fprintf('\n \t TOL      : %1.2e',atol);
    fprintf('\n \t ..... adaptive .....');
    fprintf('\n \t iters    : %8d',adapt_iters); 
    fprintf('\n \t triangles: %8d',length(adapt_L1_vertices));
    fprintf('\n \t cpu time : %1.2e',cpu_time); 
    fprintf ("\n \t .......... errors........... ") ;
    fprintf ("\n \t I  : %1.5e ", I) ;
    fprintf ("\n \t Ia : %1.5e ", Iadapt) ;
    fprintf ("\n \t RE : %1.5e \n ", RE) ;
    fprintf ("\n \t -----------------------------\n") ;
end



function plot_errors(ii,n,reV,log_re)

if ii <= 5
    switch ii
        case 1
            plotstr='m+'; linestr='m-';
        case 2
            plotstr='g+'; linestr='g-';
        case 3
            plotstr='r+'; linestr='r-';
        case 4
            plotstr='b+'; linestr='b-';
        case 5
            plotstr='c+'; linestr='b-';
        otherwise
            plotstr='k.'; linestr='w.';
    end
    %     semilogy(1:number_experiments,reV,plotstr); hold on;
    %     semilogy(1:number_experiments,log_reV(ii)*ones(1,number_experiments),...
    %         linestr,'Linewidth',3);
    semilogy(n*ones(size(reV)),reV,plotstr,'LineWidth',2); hold on;
    semilogy(n,log_re, 'ko','MarkerSize',30,'MarkerEdgeColor','k',...
        'LineWidth',2);
else
    semilogy(n*ones(size(reV)),reV,'color',rand(1,3)); hold on;
    semilogy(n,log_re, 'ko','MarkerSize',29,'MarkerEdgeColor','k');
end

