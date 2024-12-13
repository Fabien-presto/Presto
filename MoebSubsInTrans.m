function [num,den]=MoebSubsInTrans(num,den)
[imax,jmax]=size(num);
den=MoebSubs(den);
for k=1:imax
    for j=1:jmax
        num{k,j}=[zeros(1,length(den)-length(num{k,j})),num{k,j}];
        num{k,j}=MoebSubs(num{k,j});
    end
end
 
% $Id: MoebSubsInTrans.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
