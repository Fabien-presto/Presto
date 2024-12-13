function [H]=GenerateMoebMatrix(theta1,theta2,n_cf,n)

n_sample=2000;
theta=linspace(theta1,theta2,n_sample);
z=exp(i*theta);
H=[];
for m=[0:n-1]
	if (m==0)
		theta_c=linspace(theta2,2*pi+theta1,n_sample);
		zz=exp(i*theta_c);
		g=-1./MoebDen(zz);
		H=[transpose(FourierCoeffs(theta_c,g,[1:n_cf]))];
	else
		g=(i^m*MoebDen(z).^(m-1))./(MoebNum(z).^m);
		H=[H,transpose(FourierCoeffs(theta,g,[1:n_cf]))];
	end
end	



 
% $Id: GenerateDivMoebMatrix.m,v 1.2 2002/09/09 15:47:14 fseyfert Exp $ 
