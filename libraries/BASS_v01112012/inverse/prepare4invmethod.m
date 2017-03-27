function varout = prepare4invmethod(fname_meth,M,A0,figShow,opts)



if nargin<5
    opts = [];
end

try Sigma_E = opts.Sigma_E; catch Sigma_E = speye(size(M,1)); disp('Using diag noise covariance'); end
try 
    MatDir = opts.save.path;
catch
    MatDir = fullfile('Mat_files','UnknownData',fname_meth);
end

try fname_Mat = opts.save.fname_Mat; catch fname_Mat = 'Test'; end;

if isfield(opts,'Phi'), Phi = opts.Phi; else Phi = eye(1); end
if isfield(opts,'Psi'), Psi = opts.Psi; else Psi = eye(size(A0,2)); end

if strcmp(fname_meth(end-3:end),'_svd')
    fname_method = fname_meth(1:end-4);
    [U,L,V] = svd(M,'econ');
%     Phi = V(:,1:Nk)';
    Phi = V';
else
    fname_method = fname_meth;
end

try fig_title = opts.fig.title; catch fig_title = ' '; end

if ~isfield(opts,'maxIter'), opts.maxIter = 100; end
if isfield(opts,'flag_calcMetrics'), flag_calcMetrics = opts.flag_calcMetrics; else flag_calcMetrics = true; end

[Nc Nd] = size(A0);
Nk = size(Phi,1);
Ns = size(M,2);

opts_inv = opts;


%% Call specific method of choice
switch fname_method
    
    case {'reg_bf'}     %Conventionally regularized beamforming
        try temporal_basis = opts.temporal_basis; catch temporal_basis='Not'; end
        
        scaleM = 1e6;
        scaleA = 1e6;
        M = M*scaleM;
        A0 = A0*scaleA;

        tic
        reg = 1e-8;
        [S,W]=reg_bf(M,A0,reg);        
        t =toc;
        
        %Output
        varout = struct('S',S,'W',W,'time',t);
        varout.scaleM = scaleM;
        varout.scaleA = scaleA;
        
    case {'reg_bf_rank'}     %Conventionally regularized beamforming
        try temporal_basis = opts.temporal_basis; catch temporal_basis='Not'; end

        scaleM = 1e6;
        scaleA = 1e6;
        M = M*scaleM;
        A0 = A0*scaleA;

        ref_channels = 15:16;
        A = A0;
        Arec=A-repmat(mean(A(ref_channels,:),1),[size(A,1) 1]);
        Arec(ref_channels,:) = [];
        A0=Arec;
        
        M = bsxfun(@minus,M,mean(M(ref_channels,:),1) );
        M(ref_channels,:) = [];
        
        tic
        reg = 1e-8;
        [S,W]=reg_bf(M,A0,reg);        
        t =toc;
        
        %Output
        varout = struct('S',S,'W',W,'time',t);    
        varout.scaleM = scaleM;
        varout.scaleA = scaleA;
        
    
    case {'BayesMN'}
        tic
            opts_inv.Sigma_e = Sigma_E;
            opts_inv.flag_update.S = true;
            opts_inv.flag_update.alpha = true;
            opts_inv.flag_update.beta = true;
            opts_inv.flag_update.hyp_trace = true; %NB only true for MN

            [S W invAlpha invBeta,opts_out] = BayesMNandLORETA_svd(M,A0,opts_inv);            
        t =toc;
        
        %Output
        varout = struct('S',S,'invAlpha',invAlpha,'invBeta',invBeta,'time',t);
        varout.ite = opts_inv.maxIter;
        if isfield(opts_out,'logEv'), varout.logEv = opts_out.logEv; else varout.logEv = NaN; end
        
    case {'BayesLORETA','BayesLORETAtest'}
        tic
            opts_inv.maxIter = opts.maxIter;
            opts_inv.Sigma_e = Sigma_E;
            opts_inv.K = opts.K;
            opts_inv.invK = opts.invK;
            [S W invAlpha invBeta,opts_out] = BayesMNandLORETA_svd(M,A0,opts_inv);            
        t =toc;
        
        %Output
        varout = struct('S',S,'invAlpha',invAlpha,'invBeta',invBeta,'time',t);
        varout.W = W;
        varout.ite = opts_inv.maxIter;
        if isfield(opts_out,'logEv'), varout.logEv = opts_out.logEv; else varout.logEv = NaN; end

        
    case {'champagne_2'}
        try temporal_basis = opts.temporal_basis; catch temporal_basis='Not'; end
        
        tic
            flag1 = 1; flag2 = 0; flag3 = 0;
            flag4 = 0; %Do not calculate log-evidence
            [mu,dmu,ite,gamma,logEv] = feval(str2func(fname_method),...
                    A0,M,full(Sigma_E),opts.maxIter,flag1,flag2,flag3,flag4);
            S = reshape(mu,Nd,Ns);
        t =toc;
        
        %Output
        varout = struct('S',S,'gamma',gamma,'time',t);
        varout.ite = ite;
        varout.logEv = logEv;
        figure, plot(1:ite,logEv),grid on, title('Champagne')
        
        
    case {'awsm_champ'}
        tic
        if isfield(opts,'scaleMA')
            scaleM = 1e6;
            scaleA = 1e6;
            M = M*scaleM;
            A0 = A0*scaleA;
        else
            scaleM = 1;
            scaleA = 1;
        end
        maxIter = opts.maxIter;
        [gamma,S,w,sigu,like]=awsm_champ(M,A0,Sigma_E,maxIter,1,0,0,[],1);
%         [gamma,S,w,sigu,like]=awsm_champ(M,A0,[],opts.maxIter,1,0,1,[]);
        t=toc;
        
        varout = struct('S',S,'gamma',gamma,'time',t);
        varout.logEv = like;
        varout.Sigma_E = sigu;
        varout.w = w;
        varout.scaleM = scaleM;
        varout.scaleA = scaleA;
        
    case {'AquavitSpace'}
        
%         flag_update = 2;    %MacKay
        %Default Wipf-update
        if isfield(opts,'flag_update'), flag_update = opts.flag_update; else flag_update = 3; end

        tic
        [subS W gamma alpha ite] = aquavit_v3(M,A0*Psi,Phi,Sigma_E,opts);
        subSw = W'*Phi;
        subV = subS-subSw;
        S = Psi * subS;
        Sw = Psi * subSw;
        V = S - Sw;

        t=toc;

        varout = struct('S',S,'gamma',gamma,'W',W,'alpha',alpha,'time',t);
        varout.V = V;
        varout.Sw = Sw;
        varout.Psi = Psi;
        varout.Phi = Phi;
        varout.ite = ite;
        
        
        
    case {'SBLwSTE','SBLwSTEfactwise','sAquavit'}
        tic
        [V,W,invgamma,invalpha,ite] = feval(str2func(fname_method),M,A0,Psi,Phi,Sigma_E,opts);
        S = Psi*W*Phi + V;
        t=toc;
        
        varout = struct('V',V,'invgamma',invgamma,'W',W,'invalpha',invalpha,'time',t);
        varout.S = S;
        varout.Psi = Psi;
        varout.Phi = Phi;
        varout.ite = ite;
%         varout.logEv = like;
%         varout.Sigma_E = sigu;
%         varout.scaleM = scaleM;
%         varout.scaleA = scaleA;

        




%     case {'awsm_champ_rank'}
% %         A0(end-1:end,:) = [];
%         scaleM = 1e6;
%         scaleA = 1e6;
%         M = M*scaleM;
%         A0 = A0*scaleA;
% 
%         ref_channels = 15:16;
%         A = A0;
%         Arec=A-repmat(mean(A(ref_channels,:),1),[size(A,1) 1]);
%         Arec(ref_channels,:) = [];
%         A0=Arec;
%         
%         M = bsxfun(@minus,M,mean(M(ref_channels,:),1) );
%         M(ref_channels,:) = [];
%         
%         tic
%         [gamma,S,w,sigu,like]=awsm_champ(M,A0,[],opts.maxIter,1,0,1,[]);
%         t=toc;
%         
%         varout = struct('S',S,'gamma',gamma,'time',t);
%         varout.logEv = like;
%         varout.Sigma_E = sigu;
%         varout.w = w;
%         varout.scaleM = scaleM;
%         varout.scaleA = scaleA;
        
                
    case {'VB_GLM_spatio_temporal','VB_GLM_spatio_temporal_fast','Barreto_STortho1','Barreto_STortho2'}
        tic
%             [S,W,alpha,lambda,sigma,b,c] = feval(str2func(fname_method),...
%                 M,A0,Phi',opts);
            [S,W,alpha,lambda,sigma,b,c] = feval(str2func(fname_method),...
                Sigma_E^(-1/2)*M,Sigma_E^(-1/2)*A0,Phi',opts);
        t=toc;
        varout = struct('S',S,'G',W,'alpha',alpha,'time',t,...
            'lambda',lambda,'sigma',sigma,'b',b,'c',c);
        
    case {'Bolstad_ESP2_TBsparse','Bolstad_ESP2_TBfull'}
        tic
        temporal.T = Phi';
        spatial.S = Psi;
        if strcmp(fname_method,'Bolstad_ESP2_TBsparse')
            temporal.J = size(Phi,1);
            spatial.I = size(Psi,2);
        else
            temporal.J = 1;
            spatial.I = 1;
        end
        
        opts_esp = struct('maxIter',opts.maxIter,'temporal',temporal,'spatial',spatial);
%               .lambda: default according to eq.(15) in [1].
        [G S lambda CardA pwiseApprox lambda_arr lmax GAll] = Bolstad_ESP2(Sigma_E^(-1/2)*M,Sigma_E^(-1/2)*A0,opts_esp);
        t=toc;
        
        %Output
        varout = struct('S',S,'G',G,'time',t);
%         varout.ite = length(rho)-1;
        varout.lambda = lambda;
        varout.lambda_arr = lambda_arr;
        varout.lmax = lmax;
        varout.CardA = CardA;
        varout.pwiseApprox = pwiseApprox;
        varout.GAll = GAll;
        
    case {'sLORETA'}
        [Nc Nd] = size(A0);
        S = nut_sLORETA(reshape(A0,[Nc 1 Nd]),struct('y',M));
        varout.S = S;
        
    case {'True','true'}
        try S = opts.sim.S0; catch disp('True sources not specified in opts.sim.'), return; end
        varout.S = S;
        
    otherwise
        try S = opts.S; catch disp('Estimated sources needs to be specified in opts.S'), return; end
end

varout.Phi = Phi;


%% Save variables to mat-file
if ~isdir(MatDir)
    mkdir(MatDir)
end

try Sref = opts.sim.S0; catch Sref = NaN; end;
try NdFull = opts.sim.NdFull; catch NdFull = size(S,1); end

if figShow
    MatDirFig = fullfile(MatDir,'Figures');
    if ~isdir(MatDirFig)
        mkdir(MatDirFig)
    end

    figure, imagesc(S), colorbar, title(sprintf('Sources: %s',fname_meth))
    figname = fullfile(MatDirFig,[fname_Mat '_imagesc.fig']);
    saveas(gcf,figname,'fig')
    close

    try icols = opts.icols; catch icols = 1:size(S,1); end;
    clear opts
    
%     keyboard
    
    [Fgraph vert face] = glass_brain(M,A0,S,icols,NdFull,fname_meth);
    vertOr = vert;
    
    vert = vert(icols,:);
    varout.add.vert = vert;
    title(fig_title)
    figname = fullfile(MatDirFig,[fname_Mat '_glass_brain.fig']);
    saveas(gcf,figname,'fig')
    
    fig_name = fullfile('Figures/SinesSVD3',[fname_method '_glas_brain']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc',[fig_name '.eps'])
%     close


flag_fig_doc=1;

if flag_fig_doc

%     keyboard
    if size(S,1)<size(vertOr,1)  %Tilfoej en raekke med nuller pga
                                 %oprindelige fejl i subsampling hvor offset
                                 %var sat 1 i stedet for 0.
        Splot = [zeros(1,size(S,2)); S];    %Det er plads nr.1.
    else
        Splot = S;
    end
    
    fig_path_doc = 'Figures/SinesMix_MLMI_final';
    
    Svar = var(Splot,[],2);    
    [handles iDip DipColor] = plot_3DbrainAndEllipsoidal(vertOr,face,[],Svar);
    
    set(gcf,'color','white')    
%     fig_path_doc = 'c:\cs\My Pictures\BrainImages\MLMI2010_figEdit';
%     fname_method = ['ColorEdit_' fname_method];
    
    view(-90,0)
    fig_name = fullfile(fig_path_doc,[fname_method '_3DbrainAddSpheres_view1']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc',[fig_name '.eps'])
    
    view(-90,90)
    fig_name = fullfile(fig_path_doc,[fname_method '_3DbrainAddSpheres_view2']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

%     view(1,-90)
    view(-90,-90)
    fig_name = fullfile(fig_path_doc,[fname_method '_3DbrainAddSpheres_view3']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

    view(90,0)
    fig_name = fullfile(fig_path_doc,[fname_method '_3DbrainAddSpheres_view4']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc',[fig_name '.eps'])
    
    view(90,90)
    fig_name = fullfile(fig_path_doc,[fname_method '_3DbrainAddSpheres_view5']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

    view(90,-90)
    fig_name = fullfile(fig_path_doc,[fname_method '_3DbrainAddSpheres_view6']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

    
    % Show associated time series
    Sdip = Splot(iDip,:);
    Ndip = length(iDip);
    fs = 1e3;
    taxis = (0:size(S,2)-1)/fs;       %time axis

%     figure
%     hold on
%     for i=1:Ndip
%         ht(i) = plot(taxis,Sdip(i,:),'Color',DipColor(i,:));
%     end
%     hold off
%     grid on
%     xlabel('Time [s]')
    
    figure
    for i=1:Ndip
        subplot(Ndip,1,i)
        ht(i) = plot(taxis,Sdip(i,:),'Color',DipColor(i,:));
        grid on
        if i<Ndip
            set(gca,'Xtick',[]);
        end
        legend(['#' num2str(i)])
    end
    xlabel('Time [s]')
    fig_name = fullfile(fig_path_doc,[fname_method '_TimeSeriesAll']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])


    for i=1:Ndip
        figure
        ht(i) = plot(taxis,Sdip(i,:),'Color',DipColor(i,:));
        grid on
        if i<Ndip
            set(gca,'Xtick',[]);
        else
            xlabel('Time [s]')
        end
        legend(['#' num2str(i)])
        fig_name = fullfile(fig_path_doc,[fname_method '_TimeSeries_iDip' num2str(i)]);
        saveas(gcf,[fig_name '.fig'])
        print(gcf,'-depsc2',[fig_name '.eps'])
    end    
%     figure
%     hold on
%     offset = 20*Ndip;
%     for i=1:Ndip
%         offset = offset - 20;
%         ht(i) = plot(taxis,Sdip(i,:)+offset,'Color',DipColor(i,:));
%     end
%     hold off
%     grid on
%     xlabel('Time [s]')

    %Plot 3D brain at a time-shot
    tplot = 159
    [handles] = plot_3Dbrain(vertOr,face,Splot(:,tplot),struct('cMAP',jet(256),'crange',[-10 10]) );
    view(-90,0)
    fig_name = fullfile(fig_path_doc,[fname_method '_3Dbrain_view1']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

    view(-90,90)
    fig_name = fullfile(fig_path_doc,[fname_method '_3Dbrain_view2']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

%     view(1,-90)
    view(-90,-89)
    fig_name = fullfile(fig_path_doc,[fname_method '_3Dbrain_view3']);
    saveas(gcf,[fig_name '.fig'])
    print(gcf,'-depsc2',[fig_name '.eps'])

end   
 
end
% keyboard

if flag_calcMetrics
    opts_metrics.tselect = 1:size(S,2);
    try
        opts_metrics.vert = vert;
        list_active = find(sum(Sref.^2,2));
        if ~isempty(list_active), opts_metrics.vert_ref_center = vert(list_active,:); end    
    catch
    %     disp('Vertices not extracted');
    end
    Res = calcMetrics_UCSF(M,A0,S,Sref,opts_metrics);
    varout.metrics = Res;
end
    
fname_Mat = fullfile(MatDir,[fname_Mat '.mat']);
save(fname_Mat)

%EOF
end