function [val,g,H]=OptCausalLagrangian(lmult,Ls)
% Only a wrapper for EvaluateCausalLagrangian used for minimsation purposes

[val,v,g,H]=EvaluateCausalLagrangian(lmult,Ls,0);
val=-val;
g=-g;
H=-H;

% $Id: OptCausalLagrangian.m,v 1.1 2003/10/01 16:48:48 fseyfert Exp $