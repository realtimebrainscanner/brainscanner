function [Aprime HR FR hits false_alarm] = calcAprime(S,pos,pos0,opts)
%==========================================================================
% Filename: calcAprime.m (function).
% 
% Description:  Calculates the A' metric according to [1]
%
% Usage:        CC = calcCorrCoef(Strue,S,opts)
%
% Input:        S: Estimated dipole activities, Nd x Nt (dipole x time)
%               pos: 3D positions of sources in S (Nd x 3) - in [mm]
%               pos0: 3D positions of cluster centers (Ncluster x 3) - in [mm]
%               opts
%                   .thresh:
%                   .dist: Distance
%
% Output:       Aprime: HR/FR
%               HR: Hit rate
%               FR: False positive rate
%
% History:
%   - Created:  09/03/2010
%   - Modified:  
%
% Special remarks:
%       NB beregningen af Nfalse er sat til 5% af Nd - dette er dog
%       afhaengig af oploesning og vil paavirke Aprime maalet.
%
%    [1]: Snodgrass, J. & Corwin, J. "Pragmatics of measuring recognition
%         memory: Applications to dementia and amnesia" Journal of
%         Experimental Psychology: General, 1988, 117, 34-50.
%
% Copyright (C) Carsten Stahlhut, DTU IMM 2010
%==========================================================================

Ncluster = size(pos0,1);
Nd = size(S,1);

try thresh = opts.thresh; catch thresh = 0.01; end
try Dneigh = opts.dist_neigh; catch Dneigh=20; end
try D = opts.dist; catch D=10; end

Dcluster = zeros(1,size(pos0,1));
for i=1:size(pos0,1)
    Edist = calcDist(pos0((1:size(pos0,1)~=i),:),pos0(i,:));
    Dcluster(i) = min(Edist);
end

P = sum(S.^2,2);
local_peak = false(Nd,1);
hits = zeros(1,Ncluster);
false_alarm = 0;

Dneigh = min(Dneigh,min(Dcluster));

for i=1:Nd
    Edist = calcDist(pos,pos(i,:));
    ineigh = find(Edist<Dneigh);
    Neighbor(i).list = ineigh;
    Neighbor(i).Edist = Edist;

    %Find local peaks
    [dummy imax] = max(P(ineigh));
    if (ineigh(imax)==i && max(abs(S(i,:)))>max(abs(S(:)))*thresh)
        local_peak(i) = true;
        
        %Check for Hit
        Edist = calcDist(pos0,pos(i,:));
        [minDist icluster] = min(Edist);
        if D>minDist
            hits(icluster) = hits(icluster)+1;
        else
            false_alarm = false_alarm+1;
        end
    end
    
end

Nfalse = Nd*0.05;

HR = sum(hits)/Ncluster;
FR = false_alarm/Nfalse;

if FR>HR
    Aprime = 0.5 - ((FR-HR)*(1+FR-HR))/(4*FR*(1-HR));
else
    Aprime = 0.5 + ((HR-FR)*(1+HR-FR))/(4*HR*(1-FR));
end

end


function Edist = calcDist(XYZ,iDipAc)
%==========================================================================
% Filename: Edist2ActiveDip.m (function).
% 
% Description:  
%
% Input:        XYZ:    Coordinates in 3D-space
%               iADip:  Vector with indices
%
% Output:       Edist:  Euclidean distance between dipoles.
%
% History:
%   - Created:  29/05/2008
%   - Modified:  
%
% Special remarks:
%
% Copyright (C) Carsten Stahlhut (s022205), IMM DTU 2008
%==========================================================================
Nd = size(XYZ,1);
[Nrow Ncol] = size(iDipAc);
if ((Nrow>1) || (Ncol>1))
    XYZ2 = iDipAc;      %iDipAc keeps 
    NdAc = Nrow;
    Edist = inf*ones(Nd,NdAc);
    vOnes = ones(Nd,1);
    for i=1:NdAc
        Edist(:,i) = sqrt(sum((XYZ - vOnes*XYZ2(i,:)).^2,2));
    end
else
    NdAc = length(iDipAc);
    Edist = inf*ones(Nd,NdAc);
    vOnes = ones(Nd,1);
    for i=1:length(iDipAc)
        Edist(:,i) = sqrt(sum((XYZ - vOnes*XYZ(iDipAc(i),:)).^2,2));
    end
end

end
