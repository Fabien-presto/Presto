function [f]=FourierCoeffs(theta,y,k_range)
theta=reshape(theta,1,length(theta));
y=reshape(y,1,length(y));
l=length(theta);
eps0=(theta(2:l)-theta(1:l-1))/2;
eps1= [eps0,0];
eps2= [0,eps0];
eps=(eps1+eps2);
yw =eps.*y;
f=[];
v0=exp(-k_range(1)*i*theta);
vd=exp(-i*theta);
for k=k_range
    %c=sum(exp(-k*i*theta).*yw)/2/pi;
    c=sum(v0.*yw)/2/pi;
    f=[f,c];
    v0=v0.*vd;
end
 
% $Id: FourierCoeffs.m,v 1.3 2002/08/30 16:05:00 fseyfert Exp $ 
