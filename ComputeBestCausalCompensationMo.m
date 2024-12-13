function [delay,svdmin]=ComputeBestCausalCompensation(w,y,d_range,poly_order,deg_rational,n_cf,w_lim,modulus_factor)
% Computes best delay compensation by looking for best causal completion 
% delay is the optimal delay found and ecmin is a measure of the non causal part

number_of_spline_points=300;
ws=linspace(min(w),max(w),number_of_spline_points);
ys=spline(w,y,ws);
%ws=w;
%ys=y;
k=0;
e=zeros(size(d_range));
for r=d_range
    k=k+1;
    yc=ys.*exp(i*r*ws);
    [p,q,e(k)]=l2ratappeig(ws,yc,deg_rational,deg_rational);
end
clf;
plot(d_range,e);
pause;
ipos=find(d_range>0);
[svdmin,i_min]=min(e(ipos))
d_range_p=d_range(ipos);
delay=d_range_p(i_min);
