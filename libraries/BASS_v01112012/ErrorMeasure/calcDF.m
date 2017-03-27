function DF = calcDF(S,Sref,iSim,t)
%Calculation of degree of focalization
%Ref Grova et al. (2006) and Im et al. (2003)
%
% iSim should be a vector indicating with numbers indicating which
% simulated region it corresponds.
% I.e. 0: if not active
%      1: region number one
%      2: region number two
%      etc

%DF is a vector with the degrees of focalization for each of the simulated
%regions

Nreg = max(iSim);
DF = NaN*ones(Nreg,1);
for ireg = 1:Nreg
    iactive = iSim==ireg;
    
    DF(ireg) = sum( (S(iactive,t)-Sref(iactive,t)).^2)/...
         sum(Sref(iactive,t).^2);
end