function S=FreqResp(p,q,w)
s=size(p);
S=zeros([s,length(w)]);
w=i*(real(w)+imag(w));
ev_q=polyval(q,w);
for k=1:s(1)
    for j=1:s(2)
        S(k,j,:)=polyval(p{k,j},w)./ev_q;
    end
end
 
% $Id: FreqResp.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
