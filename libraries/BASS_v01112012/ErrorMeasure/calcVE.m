function VE = calcVE(M,Mrecon)
% Calculaiton of variance accounted for (%)

SSR = sum(var((M - Mrecon),0,2));
SST = sum(var(M,0,2));

% SSR = sum((M(:)-Mrecon(:)).^2);
% SST = sum(M(:).^2);
 
% assess accuracy; signal to noise (over sources)
%======================================================================
VE = 100*(SST - SSR)/SST;   %variance accounted for (%)