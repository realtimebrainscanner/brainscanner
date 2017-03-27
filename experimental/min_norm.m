function [M] = min_norm(beta, Phi, T)

L=ones(size(Phi));
alpha = 0.1;
PhiInv = (Phi*Phi' + alpha*L*L')\(Phi'*L');

x = PhiInv * T;

end