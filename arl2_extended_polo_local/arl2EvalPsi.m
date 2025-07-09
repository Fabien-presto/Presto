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
%line 211 <arl2.nw>
function [psi, Dpsi, D2psi, sys] = arl2EvalPsi(X, chart, inputData, force)

if nargin<4, force = 0; end

needGrad = (nargout>1);

Dpsi  = [];
D2psi = [];

% Expected parameters and additionnals
[p,n] = size(chart.u);
N = 2*n*p;
if length(X)<N
 error('Internal error: chart and parameters sizes mismatch')
else
 Xad = X(N+1 : end);
 X   = X(1:N);
end

% Special case: data need rotation (ALCATEL FILTERS)
switch length(Xad)
case 0
  needRotate = 0;
  data = inputData;
case 1
  needRotate = 1;
  data = feval(inputData.addProcess, inputData, Xad(1) );
otherwise
  errror('Internal error: unexpected additionnal parameters')
end

CSP.chart = chart;
[CSP.y, CSP.D] = arl2X2SP(X,n,p);

if ~force & ~arl2IsCSP(CSP)
 psi     = Inf;
 return
end

if needGrad
 [allpass,Slist,Ulist,Vlist] = arl2CSP2Allpass(CSP);
 [J, sys, dJ] = arl2UserFunction(allpass, data);
 psi = J;
 dy = arl2CSP2GradSP(CSP,dJ,Slist,Ulist,Vlist);
 dy = dy(:); Dpsi = [ real(dy) ; imag(dy) ];
 if needRotate
   [D, G] =  feval(inputData.addProcess, inputData, Xad(1) , allpass);
   Dpsi = [Dpsi ; G(:) ]; 
 end
% arl2CheckDerivative(chart, dy, X, data);
else
 allpass = arl2CSP2Allpass(CSP);
 [psi , sys]   = arl2UserFunction(allpass, data);
end
%line 266 <arl2.nw>
function arl2CheckDerivative(chart, dy, X, data)
 [p,n] = size(chart.u);
 CSP.chart = chart;
 GRAD = zeros(size(X));

%% FINITE DIFFERENCES
 DX = 1.e-5;
 % J0
 [CSP.y, CSP.D] = arl2X2SP(X,n,p);
 allpass = arl2CSP2Allpass(CSP);
 J0 = arl2UserFunction(allpass, data);

 % Perturbation k of X -> Jk
 for k = 1:length(X)
  Xk = X;
  Xk(k) = Xk(k) + DX;
  [CSP.y, CSP.D] = arl2X2SP(Xk,n,p);
  allpass = arl2CSP2Allpass(CSP);
  Jk = arl2UserFunction(allpass, data);
  GRAD(k) = (Jk - J0)/DX;
 end

%% EXPLICIT GRADIENT
 Y = dy(:);
 grad = [ real(Y) ; imag(Y)];

%% COMPARISON
 fprintf('arl2CheckDerivative: %f\n',norm(GRAD-grad));
T = 1:length(GRAD);
clf;  plot(T,GRAD,'b',T,grad,'r'); drawnow
