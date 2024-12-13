function [S]=RatApp(S,n); 
addpath /net/miaou/lib/arl2/Matlab ;
addpath ../arl2polo;
assert(isfield(S,'fourier'),'Compute fourier coeffs first');
for k=1:2
    for j=1:2
        CoeffPlus(k,j,:)=conj(S.fourier{k,j}(2:end));        
    end
end
data.value = CoeffPlus;
data.type  = 'coef';
sys0 = arl2CreateStableSystem(n,2,2);
options = { 			...
  'TolFun',		eps,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'on'   ...
};
% Schur parameters monitoring 
arl2InitMonitor([], 10, 40, 4); 
[sys, report] = arl2(10, data, sys0, options);
close all;
S.rat_sys=sys;
%resp = SSresp(sys,1000);
%theta=linspace(0,2*pi,100);
%v=conj((freqresp(ss(sys.a,sys.b,sys.c,sys.d,1),theta)))./exp(i*transpose(theta));