function [p]=ComputeTchebyPoints(xmin,xmax,n)

% Compute tcheby points including xmin and xmax (except when n=1)
if (n==1)
    p=(xmax-xmin)/2;
else
    k=[0:n-1];
    %p=cos(pi*(2*k+1)/(2*n));
    p=cos(pi*(k)/(n-1));
    a=(xmax-xmin)/(2);
    p=xmin+a*(p+1);
    p=p(end:-1:1);
end
