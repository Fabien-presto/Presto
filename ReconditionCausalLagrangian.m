function [L,Pc]=ReconditionCausalLagrangian(L,x1,x2)
% Reconditonne le lagrangien en changeant de base pour la completion
% La nouvelle base est est constituee de poly. de Chebychev sur le segment [x1,x2]
% Pc est la matrice de changement de base

s=size(L.D);
n=s(2)-1;
Paux=transpose(ChebychevBasis(x1,x2,n));
Pc=Paux(end:-1:1,:);

L.D=L.D*Pc;
L.B=L.B*Pc;

if (isfield(L,'BB'))
    L.BB=L.B'*L.B;
end
