function [p]=MyPolyfit(x,y,n,n_prescibed_zero_in_zero)
if (n_prescibed_zero_in_zero==0)
    p=polyfit(x,y,n);
else
    assert(n>=n_prescibed_zero_in_zero,'Prescribed zero at 0 should be less or equal total degree');
    % Build matrix
    nx=length(x);
    y=reshape(y,nx,1);
    A=zeros(nx,n+1-n_prescibed_zero_in_zero);
    for k=1:nx
        A(k,n-n_prescibed_zero_in_zero+1)=x(k)^(n_prescibed_zero_in_zero);
        for l=(n-n_prescibed_zero_in_zero):-1:1
            A(k,l)=A(k,l+1)*x(k);
        end
    end
    tA=transpose(conj(A));
    p=(tA*A)\(tA*y);
    p=[transpose(p),zeros(1,n_prescibed_zero_in_zero)];
end 
% $Id: MyPolyfit.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
