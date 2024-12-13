function [f,ValuesInOne]=DivideByZminusOne(fourier)
if (nargin<2)
    n_div=1;
end
[imax,jmax]=size(fourier);
ValuesInOne=zeros(imax,jmax);
f=cell(imax,jmax);
for k=1:imax
    for j=1:jmax
        ValuesInOne(k,j)=sum(fourier{k,j});
        fourier{k,j}(1)=fourier{k,j}(1)-ValuesInOne(k,j);
        n=length(fourier{k,j});
        f{k,j}=zeros(1,n);
        for l=2:n
            f{k,j}(l)=f{k,j}(l-1)+fourier{k,j}(l-1);
        end
    end
end
 
% $Id: DivideByZminusOne.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
