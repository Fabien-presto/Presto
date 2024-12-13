function [p,q,beta,sys]=Run_Rarl2_l2(thetas,values,sys0,n,adjust,stop_treshold)

data.type  = 'sample';

% The data are splined in 500 points - this is a magic constant for the moment
data.theta = L2C(linspace(C2L(min(thetas)),C2L(max(thetas)),500));
[imax,jmax,aux]=size(values);
s_values=zeros(imax,jmax,length(data.theta));
for j=1:imax
    for k=1:jmax
        s_values(j,k,:)=spline(thetas,squeeze(values(j,k,:)),data.theta);
    end
end
data.value=s_values;

options = { 			...
  'TolFun',	 stop_treshold,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'off'   ...
  'Display','iter', ...
};

if (adjust)
    data.moreParams = struct(...
    'X0', [0], ...
    'LB', [-Inf], ...
    'UB', [Inf] ...
    );
    data.addProcess = @dataRotate;
end

arl2CloseMonitor;
[sys, report] = s_arl2(n, data, sys0, options);
if (adjust)
    beta=report.Xad;
else
    beta=0;
end
% Build transfer from state space description
[num,den]=MyTf(sys.a,sys.b,sys.c,sys.d);

p=num;
q=den{1,1};


% $Id: 