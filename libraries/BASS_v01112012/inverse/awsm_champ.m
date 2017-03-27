

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Learn voxel covariances and noise level in Champagne
%
% 3/28/12
% 
% © 2012 Convex Imaging
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output: 
% gamma(nd,nd,nv) = voxel covariance matrices
% x(nv*nd,nt) = voxel current density 
% w(nv*nd,nk) = reconstruction filter
% sigu(nk,nk) = noise covariance matrix
% like(nem,1) = likelihood
%
% Input:
% y(nk,nt) = sensor data 
% f(nk,nv*nd) = lead field matrix
% sigu(nk,nk) = noise covariance matrix, [] for automatic initialization
% nem = maximum # of iterations
% nd = # of orientations
% vcs = voxel covariance structure: 0 = scalar, 1 = diagonal, 2 = general 
% nupd = 2: learn diagonal noise covariance, = 1: learn scalar, = 0: use input noise cov
% gamma0(nd,nd,nv) = initial voxel covariances, [] for automatic initialization 
% retx = 1: return voxel current density, = 0: don't 

%
% #s:
% nk = # of sensors
% nv = # of voxels
% nt = # of data time points



function [gamma,x,w,sigu,like]=awsm_champ(y,f,sigu,nem,nd,vcs,nupd,gamma0,retx);

if vcs==2 && nd>1
    [gamma x w sigu like]=awsm_champ_mat(y,f,sigu,nem,nd,nupd,gamma0,retx);
else
    [gamma x w sigu like]=awsm_champ_vec(y,f,sigu,nem,nd,vcs,nupd,gamma0,retx);
end

return




function [gamma,x,w,sigu,like]=awsm_champ_vec(y,f,sigu,nem,nd,vcs,nupd,gamma0,retx);

eps1=1e-8;
eps1z=1e-8;
[nk nvd]=size(f);
nv=nvd/nd;
nt=size(y,2);
nb=1000;

cyy=y*y'/nt;

% Initialize voxel variances

if isempty(gamma0)
    f2=sum(f.^2,1);
    invf2=zeros(1,nvd);
    ff=find(f2>0);
    invf2(ff)=1./f2(ff);
    w=spdiags(invf2',0,nvd,nvd)*f';

    if nt>nb
        inu0=zeros(nvd,1);
        for it=1:nb:nt-nb+1
            inu0=inu0+sum((w*y(:,it:it+nb-1)).^2,2);
        end
        inu0=(inu0+sum((w*y(:,it+nb:nt)).^2,2))/nt;
    else
        inu0=mean((w*y).^2,2);
    end
%    disp(max(inu0-mean((w*y).^2,2)));
    inu0v=mean(reshape(inu0,nd,nv),1);
    inu1=reshape(ones(nd,1)*inu0v,nvd,1);
    vvec=inu1;
else
    vvec=zeros(nvd,1);
    for iv=1:nv
        jv=(iv-1)*nd+1:iv*nd;
        if vcs==0
            vvec(jv)=mean(diag(gamma0(:,:,iv)));
        else
            vvec(jv)=diag(gamma0(:,:,iv));
        end
    end
end

% initialize noise covariance

if nupd>0
    if isempty(sigu)
        ndr=.1;
        sigu=ndr*mean(mean(y.^2))*eye(nk);
    else
        sigu=mean(diag(sigu))*eye(nk);
    end
end

% Learn voxel variances

v=zeros(nv,1);
figure;

like=zeros(nem,1);
for iem=1:nem
    vmat=spdiags(vvec,0,nvd,nvd);
    c=f*vmat*f'+sigu;
    [p d]=svd(c);
    d=max(real(diag(d)),0);
    invd=zeros(nk,1);
    ff=find(d>=eps1);
    invd(ff)=1./d(ff);
%    invd=1./d;
    invc=p*spdiags(invd,0,nk,nk)*p';
    
%    like(iem)=-.5*(sum(log(d))+nk*log(2*pi))-.5*sum(sum(y.*(invc*y)))/nt; 
    like(iem)=-.5*(sum(log(max(d,eps1)))+nk*log(2*pi))-.5*sum(sum(invc.*cyy));    
    subplot(2,2,1);plot((1:iem),like(1:iem));
    title(['Likelihood: ' int2str(iem) ' / ' int2str(nem)]);
    xlabel('iteration');
    set(gca(),'XLim',[0 iem]);
    
    fc=f'*invc;
    w=vmat*fc;
    x2=sum((w*cyy).*w,2);
    z=sum(fc.*f',2);

    if vcs==0
        x20=sum(reshape(x2,nd,nv),1);
        z0=sum(reshape(z,nd,nv),1);
        v0=(sqrt(max(z0,0))./max(z0,eps1z)).*sqrt(max(x20,0));
        vvec=reshape(ones(nd,1)*v0,nvd,1);
    else
        vvec=(sqrt(max(z,0))./max(z,eps1z)).*sqrt(max(x2,0));
    end
    
    if nupd>0
        fw=eye(nk)-f*w;
        sigy1=sum((fw*cyy).*fw,2);
        fgf=sum((fw*f*vmat).*f,2);
        ilam=sigy1+fgf;
        if nupd==1
            sigu=mean(ilam)*eye(nk);
        else
            sigu=diag(ilam);
        end
    end

    v=sum(reshape(vvec,nd,nv),1);
    subplot(2,2,2);plot((1:nv),v);
    title(['Voxel power: ' num2str(nv) ' / ' num2str(nv)]);
    xlabel('voxel index');
    set(gca(),'XLim',[1 nv]);
    drawnow
end

if retx==1
    x=w*y;
else
    x=[];
end

if nd==1
    gamma=reshape(vvec,1,1,nv);
else
    gamma=zeros(nd,nd,nv);
    for iv=1:nv
        gamma(:,:,iv)=diag(vvec((iv-1)*nd+1:iv*nd));
    end
end

return




function [gamma,x,w,sigu,like]=awsm_champ_mat(y,f,sigu,nem,nd,nupd,gamma0,retx);

eps1=1e-8;
eps1z=1e-8;
[nk nvd]=size(f);
nv=nvd/nd;
nt=size(y,2);
nb=1000;

cyy=y*y'/nt;

% Initialize voxel covariances

if isempty(gamma0)
    f2=sum(f.^2,1);
    invf2=zeros(1,nvd);
    ff=find(f2>0);
    invf2(ff)=1./f2(ff);
    w=spdiags(invf2',0,nvd,nvd)*f';
    
    if nt>nb
        inu0=zeros(nvd,1);
        for it=1:nb:nt-nb+1
            inu0=inu0+sum((w*y(:,it:it+nb-1)).^2,2);
        end
        inu0=(inu0+sum((w*y(:,it+nb:nt)).^2,2))/nt;
    else
        inu0=mean((w*y).^2,2);
    end
%    disp(max(inu0-mean((w*y).^2,2)));
    inu0v=mean(reshape(inu0,nd,nv),1);
    inu1=reshape(ones(nd,1)*inu0v,nvd,1);
    vmat=spdiags(inu1,0,nvd,nvd);
else
    vmat=sparse(nvd,nvd);
    for iv=1:nv
        jv=(iv-1)*nd+1:iv*nd;
        vmat(jv,jv)=gamma0(:,:,iv);
    end
end

% initialize noise covariance

if nupd>0
    if isempty(sigu)
        ndr=.1;
        sigu=ndr*mean(mean(y.^2))*eye(nk);
    else
        sigu=mean(diag(sigu))*eye(nk);
    end
end

% Learn voxel covariances

v=zeros(nv,1);
figure;

like=zeros(nem,1);
for iem=1:nem
    c=f*vmat*f'+sigu;
    [p d]=svd(c);
    d=max(real(diag(d)),0);
    invd=zeros(nk,1);
    ff=find(d>=eps1);
    invd(ff)=1./d(ff);
%    invd=1./d;
    invc=p*spdiags(invd,0,nk,nk)*p';

%    like(iem)=-.5*(sum(log(d))+nk*log(2*pi))-.5*sum(sum(y.*(invc*y)))/nt;
    like(iem)=-.5*(sum(log(max(d,eps1)))+nk*log(2*pi))-.5*sum(sum(invc.*cyy));    
    subplot(2,2,1);plot((1:iem),like(1:iem));
    title(['Likelihood: ' int2str(iem) ' / ' int2str(nem)]);
    xlabel('iteration');
    set(gca(),'XLim',[0 iem]);
    
    fc=f'*invc;
    w=vmat*fc;
%    x=w*y;
%    x2=mean(x.^2,2);
%    z=sum(fc.*f',2);

    for iv=1:nv
        jv=((iv-1)*nd+1:iv*nd);
        x2=w(jv,:)*cyy*w(jv,:)';
        z=fc(jv,:)*f(:,jv);
        
        [pz dz]=eig(z);
        dzd=max(real(diag(dz)),0);
        dz5=sqrt(max(dzd,0));
        z5=pz*diag(dz5)*pz';
        invdz5=1./sqrt(max(dzd,eps1z));
        invz5=pz*diag(invdz5)*pz';

        [px dx]=eig(z5*x2*z5);
        dx5=sqrt(max(real(diag(dx)),0));
        cx5=px*diag(dx5)*px';
        vmat(jv,jv)=invz5*cx5*invz5;
        v(iv)=sum(diag(vmat(jv,jv)));
    end

    if nupd>0
        fw=eye(nk)-f*w;
        sigy1=sum((fw*cyy).*fw,2);
        fgf=sum((fw*f*vmat).*f,2);
        ilam=sigy1+fgf;
        if nupd==1
            sigu=mean(ilam)*eye(nk);
        else
            sigu=diag(ilam);
        end
    end
    
    subplot(2,2,2);plot((1:nv),v);
    title(['Voxel power: ' num2str(nv) ' / ' num2str(nv)]);
    xlabel('voxel index');
    set(gca(),'XLim',[1 nv]);
    drawnow
end

if retx==1
    x=w*y;
else
    x=[];
end

gamma=zeros(nd,nd,nv);
for iv=1:nv
    jv=(iv-1)*nd+1:iv*nd;
    gamma(:,:,iv)=vmat(jv,jv);
end

return


