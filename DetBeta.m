function [beta_opt,err_opt]=DetBeta(S,plot_flag)
n_beta_test=1000;
if (nargin<2)
    plot_flag=0;
end       
S11=squeeze(S.value(1,1,:));
S12=squeeze(S.value(1,2,:));
S21=squeeze(S.value(2,1,:));
S22=squeeze(S.value(2,2,:));
l=length(S22);
beta_range=linspace(0,pi,n_beta_test);
err=[];
for beta=beta_range
    err=[err,norm(S11.*conj(S12)+exp(i*2*beta)*S21.*conj(S22))/sqrt(l)];
end
if (plot_flag)
    plot(beta_range,err);
end
[err_opt,ind]=min(err);
beta_opt=beta_range(ind); 
% $Id: DetBeta.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
