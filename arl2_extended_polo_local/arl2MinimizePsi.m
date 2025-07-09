% This file is part of Presto-HF, a matlab toolbox to identify a circuit from
% its response.
%
% SPDX-License-Identifier: AGPL-3.0-or-later
%
% Copyright 2025 by
%   Centre Inria de l'Université Côte d'Azur
%   2004, route des Lucioles
%   06902 Sophia Antipolis Cedex
%
% and by
%   Mines Paris - Armines
%   60, boulevard Saint-Michel
%   75006 Paris
%
% Contributors: Fabien Seyfert, Jean-Paul Marmorat, Martine Olivi
%
% Presto-HF is free software: you can redistribute it and/or modify it under
% the terms of the GNU Affero General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option) any
% later version.
%
% Presto-HF is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
% A PARTICULAR PURPOSE. See the GNU Affero General Public License for more
% details.
%
% You should have received a copy of the GNU Affero General Public License
% along with Presto-HF. If not, see <https://www.gnu.org/licenses/>.
%
%line 77 <arl2.nw>
function [newCSP, OUTPUT, HESSIAN, Xad, sys] = ...   
    arl2MinimizePsi(finalDegree, OPTIONS, CSP0, data)

syshandle = [];

[p,n] = size(CSP0.y);
L = 2*n*p;
%line 86 <arl2.nw>
% Set minimizaton problem
FUN     = 'arl2EvalPsi';
A       = [];
B       = [];
Aeq     = [];
Beq     = [];
LB      = -ones(L,1);
UB      =  ones(L,1);
NONLCON = 'arl2SPModulus';
%line 97 <arl2.nw>
% Starting point 
chart = CSP0.chart;
X     = arl2SP2X(CSP0.y, CSP0.D);
%line 104 <arl2.nw>
if isfield(data, 'moreParams')
 LB = [LB ; data.moreParams.LB(:) ];
 UB = [UB ; data.moreParams.UB(:) ];
 X  = [X ; data.moreParams.X0(:) ]; 
 L = length(X);
end
%line 112 <arl2.nw>
% First evaluation
FVAL  = arl2EvalPsi(X, chart, data);
%line 116 <arl2.nw>
% Minimization loop
iter = 0; iterMax = 20;
iterations = 0;
funcCount  = 0;
while iter < iterMax
 iter = iter + 1;
%line 124 <arl2.nw>
 % Call minimizing function
 warning off
 [X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN]  = ... 
    arl2Min(FUN, X, A, B, Aeq, Beq, LB, UB, NONLCON, OPTIONS, chart, data);
 warning on
 iterations = iterations + OUTPUT.iterations;
 funcCount  = funcCount  + OUTPUT.funcCount;
 if ~strcmp(OPTIONS.Display, 'off')
   try, delete(syshandle), end
%   global arl2ShowChartChanges
%   if ~isempty(arl2ShowChartChanges)
%     if arl2ShowChartChanges, syshandle = arl2ShowCurrentSys(X, chart, data); end
%  end
 end
 if FVAL < OPTIONS.TolFun ,  
%  disp('Absolute minimum ?'); break
 end
 switch EXITFLAG
%line 141 <arl2.nw>
 case 0
 % Change to adapted chart, show if verbose, and continue
 [X, chart] = arl2ChangeChart(X, chart);

 if ~strcmp(OPTIONS.Display, 'off')
   fprintf('Change to adapted chart - FVAL = %e / %e \n', ...
             FVAL, arl2EvalPsi(X,chart,data) )
 end
%line 151 <arl2.nw>
 case 1
% disp('Local minimum ?')
 % Check for local minimum
 currentDegree = length(chart.w);
 if currentDegree < finalDegree
   if ~strcmp(OPTIONS.Display, 'off')
     fprintf('Change degree - FVAL = %e / %e \n', ...
             FVAL, arl2EvalPsi(X,chart,data) );
   end
   [X , chart] = arl2IncrementDegree(X, chart);
%X
%chart.u
%chart.w
%Z = arl2X2SP(X,n+1,p)
%sum(Z.*conj(Z))
%  [psi, dpsi] = arl2EvalPsi(X,chart,data, 1) 
%length(dpsi)
%n
%p
%  Z = arl2X2SP(dpsi,n,p)
%  if ~arl2Continue, dbquit, end

else
%  disp('Maximum degree reached')
  break
 end
%line 179 <arl2.nw>
 otherwise
%line 182 <arl2.nw>
 % End switch
 end
%line 186 <arl2.nw>
% End minimization loop
end

OUTPUT.iterations    = iterations;
OUTPUT.funcCount     = funcCount;
OUTPUT.finalValue    = FVAL;
OUTPUT.TolFun        = OPTIONS.TolFun;
Xad = X(2*n*p+1:end);
%line 196 <arl2.nw>
% Build new CSP from new point X and last chart
[p,n]  = size(chart.u);
[y, D] = arl2X2SP(X,n,p);
newCSP.chart = chart;
newCSP.y = y;
newCSP.D = D;
[psi, Dpsi, D2psi, sys] = arl2EvalPsi(X, chart, data);
%line 205 <arl2.nw>
function h = arl2ShowCurrentSys(X, chart, data)
[psi, Dpsi, D2psi, sys] = arl2EvalPsi(X, chart, data);
h = SSplot(sys,'r');
