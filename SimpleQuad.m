function integral=SimpleQuad(x,y)
% Computes the integral of the function y=f(x) using the trapeze method
x=reshape(x,1,length(x));
y=reshape(y,1,length(y));
xq=([x([2:length(x)]),x(length(x))]-[x(1),x([1:length(x)-1])])/2;
integral=sum(xq.*y);
 
% $Id: SimpleQuad.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
