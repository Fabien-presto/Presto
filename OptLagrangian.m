function [val,g,H]=OptLagrangian(lmult,Ls)
% Only a wrapper for EvaluateLagrangian used for minimsation purposes

[val,g,v,ecm,ecc,ed,H]=EvaluateLagrangian(lmult,Ls);
val=-val;
g=-g;
H=-H;
 
% $Id: OptLagrangian.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
