function E = spatialAc(Q,S,plotOn)
%==========================================================================
% Filename: spatialAc.m (function).
% 
% Description:  
%
% Input:        Echoice:      SPM structure
%               Q: Structure consisting of (not all of the fields have to
%               be specified, it depends on the error measure):
%                   - EdistTable: Distance table to active dipoles.
%                     (Nd x NdActive)
%
% Output:       EdistTable: Euclidean distance table with distances to
%                           active dipoles.
%
% History:
%   - Created:  30/05/2008
%   - Modified: 17/10/2008: Changed the order of accessing EdistTable to
%                           when obtaining Edist, since EdistTable is
%                           expected to be calculated as Nd x Nd_true,
%                           where Nd is the number of dipoles in the
%                           forward model used in the reconstruction
%                           methods and Nd_true is the number of dipoles in
%                           true forward model.
%                           Note since Edist was changed the line with
%                           minEdist was also changed such that the min
%                           dist of each columns is now found.
%
% Special remarks:
%
% Copyright (C) Carsten Stahlhut (s022205), DTU IMM 2008
%==========================================================================

try,
    EdistTable = Q.EdistTable;
    iDipAc = Q.iDipAc;
catch,
    iDipAc = find(sum(abs(Q.STrue),2));
    EdistTable = Edist2ActiveDip(Q,iDipAc);
end


EnergyS = sum(S.^2,2);
EnergyS = EnergyS/max(EnergyS(:));

% iDip = find(EnergyS);
EnergyThres = 0;
iDip = (EnergyS > EnergyThres);

Edist = EdistTable(:,iDip);
minEdist(iDip) = min(Edist,[],1);
% minEdist(~iDip) = inf;

[EdistSort,iEdist] = sort(minEdist);
Nd = length(EdistSort);
Np = sum(iDip);
% cumFreq = (1:Nd)/Nd;
cumFreq = (1:Np)/Np;        %Er vist ikke korrekt!

E.dist = Edist;
E.SpatialAc.LocErr = EdistSort;
E.SpatialAc.iLocErr = iEdist;
E.SpatialAc.cumFreq = cumFreq;



% Find false negative
iDipAcLogic = false(Nd,1);
iDipAcLogic(iDipAc) = 1;

% EdistTable2 = Edist2ActiveDip(Q,find(iDip));
% Edist2 = EdistTable2(iDipAcLogic,:);
% minEdist2(iDipAcLogic) = min(Edist2,[],2);
% % minEdist(~iDip) = inf;
% 
% [Edist2Sort,iEdist2] = sort(minEdist2);
% Nd = length(Edist2Sort);
% Nn = sum(iDipAcLogic);
% % cumFreq = (1:Nd)/Nd;
% cumFreq2 = (Nn-(1:Nn))/Nn;
% 
% 
% % Nn = sum(iDipAcLogic);
% 
% % FPLE = (Np-(1:Np))/Np;
% FPLE = EdistSort/max(EdistSort);
% FPLE = (0:Np-1)/Np;
% 
% iDipDepThres = false(Nd,1);
% for jj=1:Np
%     iDipDepThres(iEdist) = 1;
%     N = find(iDip-iDipAcLogic>0);
%     FNLE(jj) = Nn-
% end


% if plotOn
%     figure
%     stairs(EdistSort,cumFreq)
%     xlabel('localization error [mm]')
%     ylabel('cumulative frequency')
%     title('spatial accuracy')
%     grid on
% end
