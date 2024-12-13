function [L]=ComputeStructureForLagrangian(x,y,xc,yc,fs,causalc,n,cflag,iso_flag)
% Builds up the structure used by EvaluateLagrangian for evaluation
% purposes
% This structure contains all the necessary parmeters to compute the
% value at a dual point lambda.

x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
xc=reshape(xc,length(xc),1);
yc=reshape(yc,length(yc),1);
L.P=GeneratePowerMatrix(x,n);
L.Pc=GeneratePowerMatrix(xc,n);
L.y=y;
L.yc=yc;
L.f=reshape(fs.f,length(fs.f),1);
% Generate two different matrices for estimation of the non-causal part
% depending of the kind of isometry used.
if (iso_flag)
	L.H=GenerateDivMoebMatrix(fs.theta1,fs.theta2,length(L.f),n+1);
else	
	L.H=GenerateMoebMatrix(fs.theta1,fs.theta2,length(L.f),n+1);
end
L.causalc=causalc;
L.norme_f2=fs.norm_f2;

% Flag concerning the type of contrained problem to be solved
L.cflag=cflag;
 
% $Id: ComputeStructureForLagrangian.m,v 1.5 2002/09/09 15:47:14 fseyfert Exp $ 
