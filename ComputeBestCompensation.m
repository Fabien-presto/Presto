function [d_opt,p,errs]=ComputeBestCompensation(w,y,d_range,poly_order,zeros_at_inf)
% Computes best delay compensation. Basically an optimized and optimized 
% least square solver for this problem

w=reshape(w,length(w),1);
y=reshape(y,length(y),1);

if (zeros_at_inf==0)
    P=GeneratePowerMatrix(1./w,poly_order);
else
    P=GeneratePowerMatrix(1./w,poly_order+zeros_at_inf);
    % Leave at least 2 freedom parameter
    i_max=max(poly_order+1,zeros_at_inf+3);
    P=P(:,[zeros_at_inf+1:i_max]);
    %P=P(:,[zeros_at_inf+1:poly_order+zeros_at_inf+1]);
end
tP=transpose(conj(P));
tPP=tP*P;
itPP=inv(tPP);
d_opt=-1;
e_opt=0;
errs=[];
for r=d_range
    yy=y.*exp(i*r*w);
    tPy=tP*(yy);
    v=itPP*tPy;
    err=norm(P*v-yy);
    errs=[errs,err];
    if (d_opt==-1 | err<e_opt)
        d_opt=r;
        e_opt=err;
        p=v;
    end
end  
p=[zeros([zeros_at_inf,1]);p];
p=transpose(p([length(p):-1:1]));

 
% $Id: ComputeBestCompensation.m,v 1.4 2003/10/01 16:40:50 fseyfert Exp $ 
