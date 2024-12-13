function [S]=ReshapeS(S,check_flag)

if (nargin < 2)
	check_flag=1;
end
n=length(S.freq);
s=size(S.value);
assert(n==s(3),'S.value should have same length as S.freq');
if (check_flag ~=-1)
	assert(min(S.freq)*max(S.freq)<0,'S.freq should have sorted frequencies around zero - S is supposed to be low pass');
end

[S.freq,ind_s]=sort(S.freq);
S.value=S.value(:,:,ind_s);
S.freq=reshape(S.freq,1,n);

 
% $Id: ReshapeS.m,v 1.5 2002/09/09 15:47:14 fseyfert Exp $ 
