function [U S V] = spsvd(X,TOL,k)
%==========================================================================
% Filename: spsvd.m (function).
% 
% Description:  Singular value decomposition with sparse representation of
%               the singular values.
%
% Input:        X:  Matrix to be decomposed
%               TOL:Relative tolerance on the singular values (inverse
%               condtion number, default is TOL=1e-6, i.e. at most
%               cond(U*S*V') = 1e6.
%
% Output:       U: Left singular vectors
%               S: Singular values (sparse matrix)
%               V: Right singular vectors
%
% History:
%   - Created:  26/02/2009
%   - Modified: 
%
% Special remarks:
%               E.g. inversion of covariance X: invX = U*(1./b)*U';
%
% Copyright (C) Carsten Stahlhut, DTU Informatics 2009
%==========================================================================

[m,n] = size(X);

if nargin<3
    k = min(m,n); 
    try dummy = TOL; catch TOL = 1e-6; end;
else
    if isempty(TOL); TOL = 1e-6; end
end

if m<=n
    decomMethod = 'econ';
else
    decomMethod = 0;
end

[U s V] = svd(X,decomMethod);
s = diag(s);

j = find(s.*[ones(k,1);zeros(length(s)-k,1)] >= max(s)*TOL);

if nargout<2
    U = s;
else
    U = U(:,j);
    V = V(:,j);
    S = spdiags(s(j),0,length(j),length(j));
end