function Phi = createTemporalBasisSet(X,method,opts)

if nargin<3
    opts = [];
end

if isfield(opts,'tol'), tol = opts.tol; else tol=1e-16; end

[Nc Ns Ntrl] = size(X);

if Ntrl>1
    X = mean(X,3);
    disp('X is 3 dimension or more. Using mean response X across trials (dim 3).')
end

switch lower(method)
    case 'svd'
        [U Sv V] = svd(X,'econ');
        sv = diag(Sv);
        sv = sv( sv > sv(1)*tol);
        Nu = length(sv);
        Phi = V(:,1:Nu)';
        
    case 'orthofilt'
        [Lo_D,Hi_D,Lo_R,Hi_R] = orthfilt(W)
    
    case 'cwt'    
        if isfield(opts,'ChNo'), ChNo = opts.ChNo; else ChNo = 1; end        
        x = X(ChNo,:);
        
        
    case 'swt'
        if isfield(opts,'wname')
            wname = opts.wname;
        else
            wname = 'sym6';
            disp('No wavelet types specified. Symlets will be used (ortogonal).')
        end
        
        if isfield(opts,'ChNo'), ChNo = opts.ChNo; else ChNo = 1; end
%         if isfield(opts,'level'), level = opts.level; else level = 8; end
        
        x = X(ChNo,:);
        level = log2(length(x));
        
        [swa, swd] = swt(x,level,wname);
        
        Phi = [swa; swd];
        
    case 'slepian'
        if isfield(opts,'fs'), fs = opts.fs; else fs = input('Specify samplingfrequency: '); end
        time_halfbandwidth = 4;
        [dps_seq,lambda] = dpss(Ns,time_halfbandwidth);
        Phi = dps_seq';
        
    case 'MP'
%         load cuspamax;
%         dict = {{'sym4',2},'sym4','dct'};
%         mpdict = wmpdictionary(length(cuspamax),'LstCpt',dict);
%         [YFIT2,R2,COEFF2,IOPT2,QUAL2] = wmpalg('OMP',cuspamax,mpdict);
        
    otherwise
        error('Not a recognized method - choose another')
end


fprintf('Number of temporal functions: %.0f\n',size(Phi,1))
