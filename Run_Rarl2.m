function [p,q,sys,beta]=Run_Rarl2(fourier_,n,adjust,stop_treshold)

% Put fourier coeffs in an array
max_size=0;
[imax,jmax]=size(fourier_);
for k=1:imax
    for j=1:jmax
        max_size=max(max_size,length(fourier_{k,j}));
    end
end
fourier=zeros(imax,jmax,max_size);
for k=1:imax
    for j=1:jmax
        l=length(fourier_{k,j});
        fourier(k,j,1:l)=fourier_{k,j};
    end
end
% Set the proper part to 0
data.valinf=zeros(imax,jmax);

data.value = fourier;
data.type  = 'coef';
% Initializes the rotation angle of S12 and S21
if (adjust)
    data.moreParams = struct(...
    'X0', [0], ...
    'LB', [-Inf], ...
    'UB', [Inf] ...
    );
    data.addProcess = @dataRotate;
end
%sys0 = arl2CreateStableSystem(n,2,2);
sys0=fc2ir(fourier,n);
options = { 			...
  'TolFun',		stop_treshold,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'off'   ...
  'Display',            'iter'
};
% Schur parameters monitoring 
%arl2InitMonitor([], 10, 40, 4); 
arl2CloseMonitor
[sys, report] = s_arl2(n, data, sys0, options);
if (adjust)
    beta=report.Xad;
else
    beta=0;
end
% Build transfer from state space description

[num,den]=MyTf(sys.a,sys.b,sys.c,sys.d);
for k=1:imax
    for j=1:jmax
        if (length(num{k,j})==n+1)
            num{k,j}(1)=[];
        end
    end
end
p=num;
q=den{1,1}; 
% $Id: Run_Rarl2.m,v 1.6 2003/12/04 15:53:04 fseyfert Exp $ 
