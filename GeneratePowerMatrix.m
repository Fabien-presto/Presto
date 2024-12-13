function P=GeneratePowerMatrix(x,n)
% Generates a matrix containing the powers of the vector x - used
% for evaluation of a polynomial at x. 

x=reshape(x,length(x),1);
P=ones(length(x),1);
for k=2:n+1
    P=[P,squeeze(P(:,k-1)).*x];
end
 
% $Id: GeneratePowerMatrix.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
