function [num,den]=MyTf(a,b,c,d)

[msg,a,b,c,d]=abcdchk(a,b,c,d); error(msg);

[n,m]=size(d);
num=cell(n,m);
den=cell(n,m);

for k=1:m
	[num_aux,den_aux]=ss2tf(a,b,c,d,k);
	for l=1:n
		num{l,k}=reshape(num_aux(l,:),1,length(den_aux));
		den{l,k}=reshape(den_aux(:),1,length(den_aux));
	end
end
