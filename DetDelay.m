function [p,delay,errmin]=DetDelay(omega_in,y_in,delay_range,poly_order,zeros_at_inf,omega_lim,opt_flag,plot_flag)

if (nargin<8)
    plot_flag=0;
end
if (nargin<7)
    opt_flag=1;
end
if (nargin<6)
    zeros_at_inf=0;
end
l=length(omega_in);
assert(l==length(y_in),'omega and y should have same length');
omega=reshape(squeeze(omega_in),l,1);
y=reshape(squeeze(y_in),l,1);
[ind]=union(find(omega>abs(omega_lim)),find(omega<-abs(omega_lim)));
omega_inf=omega(ind);
y_inf=y(ind);
err=[];
if (length(omega_inf)<2*poly_order+1)
    fprintf('Warning: approximation degree is greater than 2x number of data points \n')
    assert(length(omega_inf)>poly_order,'Approximation degree is greater than number of data points');
end
%for r=delay_range
    %p=polyfit(1./omega_inf,exp(r*i.*omega_inf).*y_inf,poly_order);
 %   p=MyPolyfit(1./omega_inf,exp(r*i.*omega_inf).*y_inf,poly_order,zeros_at_inf);
  %  err=[err,norm(polyval(p,1./omega_inf)-((exp(r*i.*omega_inf)).*y_inf))];
  %end
[delay,p,err]=ComputeBestCompensation(omega_inf,y_inf,delay_range,poly_order,zeros_at_inf); 
if (plot_flag)
    plot(delay_range,err);
end
[errmin,indm]=min(err);
if ((indm==1 | indm==length(err)) & opt_flag)
    fprintf('Warning no local minimum found in delay det. - increases delay range \n');
end
step=sqrt((max(abs(omega_inf))-abs(omega_lim))/length(omega_inf));
errmin=errmin*step;
delay=delay_range(indm);
%p=MyPolyfit(1./omega_inf,exp(delay*i.*omega_inf).*y_inf,poly_order,zeros_at_inf);
 
% $Id: DetDelay.m,v 1.3 2002/08/30 16:05:00 fseyfert Exp $ 
