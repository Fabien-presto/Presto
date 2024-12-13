function [P]=ChebychevBasis(x1,x2,n)
% Returns a Chebychev basis of chebychev polynmials of max order n on the
% segment [x1,x2]
% Stores polynomials lines by lines

P=zeros(n+1,n+1);

a=(x2-x1)/2;
b=(x2+x1)/2;

P(1,n+1)=1;
for k=1:n
    v=(pi/2+[0:1:k-1]*pi)/k;
    v=cos(v)*a+b;
    aux=poly(v);
    aux=aux/polyval(aux,a+b);
    P(k+1,[n+1-k:n+1])=aux;
end
    