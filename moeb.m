function [y]=moeb(x,inv_flag)
%y=(x+1)./(x-1);
if (nargin<2 | inv_flag==0) 
 y=MoebNum(x)./MoebDen(x);
else
 y=MoebDen(x)./MoebNum(x);
end
 
% $Id: moeb.m,v 1.4 2002/09/09 15:47:14 fseyfert Exp $ 
