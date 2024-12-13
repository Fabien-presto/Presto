function [c]=CompScalarFourierCoeffs(w,y,p_comp,f_range,n_comp,n_spline)
w=reshape(w,1,length(w));
y=reshape(y,1,length(y));
theta_i=L2C(w);
theta1=L2C(max(w));
theta2=L2C(min(w));
theta_c=linspace(theta1,2*pi,n_comp);
theta_c(1)=[];
theta_c(end)=[];
theta_cb=linspace(0,theta2,n_comp);
theta_cb(1)=[];
theta_cb(end)=[];
theta_c=[theta_c,theta_cb];
y_c=polyval(p_comp,1./C2L(theta_c));
y_t=[y,y_c];
theta_t=[theta_i,theta_c];
[theta_t,inds]=sort(theta_t);
y_t=y_t(inds);
theta_s=linspace(0,2*pi,n_spline);
y_s=spline(theta_t,y_t,theta_s);
c=FourierCoeffs(theta_s,y_s,f_range);
 
% $Id: CompScalarFourierCoeffs.m,v 1.4 2002/08/30 16:05:00 fseyfert Exp $ 
