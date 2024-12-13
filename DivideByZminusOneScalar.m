function [f]=DivideByZminusOneScalar(fourier) 
% Divides an 1/z expansion by z-1 after substraction of its value in 1 - 
%coeffs. are supposed to be sorted like f=a1+a2/z+a3/z^2 (fourier=[a1,a2,...])

f=zeros(size(fourier));
fourier(1)=fourier(1)-sum(fourier);
for l=2:length(fourier)
	f(l)=f(l-1)+fourier(l-1);
end

%$Id: DivideByZminusOneScalar.m,v 1.1 2002/09/09 15:49:47 fseyfert Exp $
