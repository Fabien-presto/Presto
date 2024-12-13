function [H]=GenerateMoebMatrix(theta1,theta2,n_cf,n)

n_sample=2000;
theta=linspace(theta1,theta2,n_sample);
z=exp(i*theta);
H=[];
for m=[0:n-1]
	if (m==0)
		v=[];
		for k=[1:n_cf]
			v=[v;(exp(-i*k*theta2)-exp(-i*k*theta1))/(-i*k*2*pi)];
		end
		H=[v];
	else
		H=[H,transpose(FourierCoeffs(theta,(i*moeb(z,1)).^m,[1:n_cf]))];
	end
end	



 
% $Id: GenerateMoebMatrix.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
