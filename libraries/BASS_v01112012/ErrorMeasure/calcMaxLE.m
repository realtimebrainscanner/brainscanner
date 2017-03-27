function [maxLE Q meanLE LEG LEGT] = calcMaxLE(Q,S,Na,vi,threshold,vert,vert_vi)
%==========================================================================
% Filename: calcMaxLE.m (function).
% 
% Description:  Calculates the maximum localization error, that is required
%               to recover at least the amount of the sources specified by
%               threshold between 0 and 1 (corresponding to 0% and 100%).
%
% Usage:        maxLE = calcMaxLE(Q,Na,threshold)
%
% Input:        Q:  Structrue with
%                       .EdistTable: Look up table with distances to the 
%                       .STrue
%                       .S
%               Na: The number of simulated dipoles, which is active.
%               vi: Indicies on dipoles the localization error (LE) should
%                   be measured with respect to. vi should keep the most
%                   active dipoles.
%               Threshold:  Level which should be recovered, [0 - 1],
%                           default is 50%.
%
% Output:       maxLE:  Maximum localization error in [mm].
%               Q:      Updated structure
%               meanLE:
%               LEG:    Localization error between gravidity of true and
%                       estimated sources.
%               LEGT:   Localization error between gravidity of true and
%                       estimated thresholded sources.
%
% History:
%   - Created:  30/05/2008
%   - Modified:  
%
% Special remarks:
%               [1]: C. Phillips, J. Mattout, M.D. Rugg, P. Maquet, and
%                    K.J. Friston, "An empirical Bayesian solution to the
%                    source reconstruction problem in EEG", NeuroImage 24,
%                    997-1011, 2005.
%
%   vi = [1936 2407];           %NB dette er de gamle placeringer!!!
%   vi = [1936 2460];           %NB dette er de nye placeringer!!!
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2008
%==========================================================================

if nargin < 5
    threshold = .50;
end

try, 
    EdistTable = Q.EdistTable;
    Q.inv = Q.inv;
catch,
    [EdistTable, D] = Edist2ActiveDip(Q,vi);
    Q.EdistTable = EdistTable;
    Q.inv = D.inv;
end

Q.iDipAc = vi;
E = spatialAc(Q,S,1);
EdistSort = E.SpatialAc.LocErr;
iEdist = E.SpatialAc.iLocErr;
Q.EdistSort = EdistSort;
Q.iEdist = iEdist;

Nest = length(iEdist);

k = find((1:Nest)/Na > threshold,1,'first');
if isempty(k)
    maxLE = inf;
else
    maxLE = EdistSort(k);
end

EnergyS = sum(S.^2,2);
[dummy ivi] = max(EnergyS);
i1maxDip = ivi;
EnergyS(ivi) = 0;
[dummy ivi] = max(EnergyS);
i2maxDip = ivi;

LEtemp = min(EdistTable([i1maxDip i2maxDip],:),[],2);
maxLE = max(LEtemp);

meanLE = sum(LEtemp)/2;


%Clustering: Find
[dist2S0, iclass] = min(EdistTable(:,vi),[],2);
clear dist2S0


%%
% Center of gravity calculated based on energy
%
Nvi = length(vi);

%Center of gravity no thresholding
Rest = NaN*one(3,Nvi);     %Gravity calculated for each time point.
LEG = Rest;
for i=length(vi)
%     Rest(i,:) = vert(repmat(iclass==i,[1 3]))'*S(iclass==i,:);
    Rest(:,i) = vert(iclass==i,:)'*EnergyS(iclass==i)/...
        sum(EnergyS(iclass==i));
    LEG(i) = sqrt((vert_vi(i,:)'-Rest(:,i)).^2);
end

% threshold = .5;
%Center of gravity based thresholding at 50 percent of maximum.
RTest = NaN*one(3,Nvi);     %Gravity calculated for each time point.
LEGT = RTest;
iclass(EnergyS/max(EnergyS) < threshold*max(EnergyS)) = 0;
for i=length(vi)
%     Rest(i,:) = vert(repmat(iclass==i,[1 3]))'*S(iclass==i,:);
    RTest(:,:,i) = vert(iclass==i,:)'*EnergyS(iclass==i)/...
        sum(EnergyS(iclass==i));
    LEGT(i) = sqrt((vert_vi(i,:)'-RTest(:,i)).^2);
end

% %%
% % Gravity at each time point
% %
% 
% %Center of gravity no thresholding
% Ns = size(S,2);
% Rest = NaN*one(3,Ns,length(vi));     %Gravity calculated for each time point.
% LEG = Rest;
% for i=length(vi)
% %     Rest(i,:) = vert(repmat(iclass==i,[1 3]))'*S(iclass==i,:);
%     Rest(:,:,i) = (vert(iclass==i,:)'*S(iclass==i,:))./...
%         repmat(sum(S(iclass==i,:),1),[3 1]);
%     LEG(:,:,i) = sqrt((repmat(vert_vi(i,:)',[1 Ns])-Rest(:,:,i)).^2);
% end
% 
% threshold = .5;
% %Center of gravity based thresholding at 50 percent of maximum.
% RTest = NaN*one(3,Ns,length(vi));     %Gravity calculated for each time point.
% LEGT = RTest;
% % iclass(EnergyS/max(EnergyS) < threshold*max(EnergyS)) = 0;
% iclass(S./max(abs(S(:))) < threshold*max(abs(S(:)))) = 0;
% for i=length(vi)
% %     Rest(i,:) = vert(repmat(iclass==i,[1 3]))'*S(iclass==i,:);
%     RTest(:,:,i) = (vert(iclass==i,:)'*S(iclass==i,:))./...
%         repmat(sum(S(iclass==i,:),1),[3 1]);
%     LEGT = sqrt((repmat(vert_vi(i,:)',[1 Ns])-RTest(:,:,i)).^2);
% end