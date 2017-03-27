

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% REG_BF: conventionally regularized beamforming 
%
% --- Version 12/28/11 -----
% 
% © 2011 Convex Imaging
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Output: 
% x(nv*nd,nt) = voxel current density 
% w(nv*nd,nk) = reconstruction filter
%
% Input:
% y(nk,nt) = sensor data (minus reference!) 
% f(nk,nv*nd) = lead field matrix
% reg = relative regularization const
% reg values in common use range between 1e-4, 1e-8, 1e-12


function [x,w]=reg_bf(y,f,reg);

[nk nvd]=size(f);
nt=size(y,2);

% Compute data covariance

cyy=y*y'/nt;
[p d]=svd(cyy);
d=diag(d);

% Regularize

ff=find(d/d(1)<reg);
if isempty(ff)
    alp=0;
else
    alp=reg*d(1);
end
%disp('alp');disp(alp/d(1));

%cyy=cyy+alp*eye(nk);
d=d+alp;

% Estimate voxel activity

invc=p*diag(1./d)*p';
w=f'*invc;
z=sum(w.*f',2);
w=((1./z)*ones(1,nk)).*w;
x=w*y;

return

