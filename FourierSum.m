function [y]=FourierSum(theta,f,first_index)
theta=reshape(theta,1,length(theta));
f=reshape(f,1,length(f));
l=length(f);
k_range=[first_index:first_index+l-1];
y=zeros(size(theta));
for k=k_range
    y=y+exp(k*i*theta)*f(k+1-first_index);
end 
% $Id: FourierSum.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
