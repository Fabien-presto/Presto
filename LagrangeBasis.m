function [P]=LagrangeBasis(x)
% Generate a basis of Lagrange polynomials according to the set of points x
% Store them line by line in a matrix p

n=length(x);
P=zeros(n,n);
for k=1:n
    aux=poly([x(1:k-1),x(k+1:n)]);
    aux=aux/polyval(aux,x(k));
    P(k,:)=aux;
end
