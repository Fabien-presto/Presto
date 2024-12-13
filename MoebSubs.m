function [q]=MoebSubs(p)
n=length(p);
q=[0];
p1=polypow([1 1],n-1);
p2=[1];
for k=1:n
     q=polyadd(q,p(k)*polymult(p1,p2));
    [p1]=deconv(p1,[1 1]);
    [p2]=conv(p2,[1 -1]);
end 
% $Id: MoebSubs.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
