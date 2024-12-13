function [H]=GenerateHankel(theta1,theta2,n,m)
% Generates the Hankel matrix associated with the linear operator
% P_{\bar H^2}(g.Chi_J) where J is the interval the subarc of the unit
% circle determined by theta1 and theta2 (trigo sens). n,m is the size
% of the resulting H     
  
H=zeros(n,m);
for k=1:n
    H(k,1)=(exp(i*k*theta2)-exp(i*k*theta1))/(k*2*pi*i);
    for j=1:min(m-1,k-1)
        H(k-j,j+1)=H(k,1);
    end
end
for l=2:m
    k=n+l-1;
    H(n,l)=(exp(i*k*theta2)-exp(i*k*theta1))/(k*2*pi*i);
    for j=l+1:m
        H(n+l-j,j)=H(n,l);
    end
end
 
% $Id: GenerateHankel.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
