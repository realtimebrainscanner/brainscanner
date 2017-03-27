function [CC list_active] = calcCorrCoef(STrue,S,opts)
%==========================================================================
% Filename: calcCorrCoef.m (function).
% 
% Description:  Calculates the correlation coefficient (CC) between the
%               true and estimated dipole activities. The CC is calculated
%               using the entire time series of all distributed dipoles, [1].
%
% Usage:        CC = calcCorrCoef(Strue,S,opts)
%
% Input:        STrue: Simulated dipole activities, Nd x Nt (dipole x time)
%               S:     Estimated dipole activities, Nd x Nt (dipole x time)
%               opts.
%                      .choice
%
% Output:       CC:    Sum of squared error score (Note its the log of the
%                      mean squared error).
%
% History:
%   - Created:  30/05/2008
%   - Modified:  
%
% Special remarks:
%               [1]: J. Daunizeau and K.J. Friston, "A mesostate-space
%                    model for EEG and MEG", NeuroImage 38, 67-81, 2007.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

try choice = opts.choice; catch choice = 'maxDip'; end

Nd = size(S,1);

if strcmp(choice,'maxDip')
    [dummy imaxT] = max(sum(STrue.^2,2));
    [dummy imax] = max(sum(S.^2,2));
    
    Cmat = corrcoef([STrue(imaxT,:)' S(imax,:)']);
%     [dummy ii] = max(abs(Cmat(:)));
    CC = Cmat(1,2);
else
    try
        list_active = opts.list_active;
    catch
        list_active = find(sum(abs(STrue),2));
    end
    P = sum(STrue(list_active,:).^2,2);
    [dummy isort] = sort(P,'descend');
    list_active = list_active(isort);
    Nactive = sum(logical(list_active));
    
    CC = nan*ones(Nactive,Nd);
%     CC = nan*ones(Nd,1);
    for iac=1:Nactive
        for ii=1:Nd
            temp = corrcoef([STrue(list_active(iac),:)' S(ii,:)']);
            CC(iac,ii) = temp(1,2);
        end
    end
end