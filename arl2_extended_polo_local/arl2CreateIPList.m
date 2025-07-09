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
%line 515 <arl2.nw>
function ipList = arl2CreateIPList(lm, maxIP, data)
%line 519 <arl2.nw>
% Check old CSP
oldCSP = lm.CSP;
oldX = arl2SP2X(oldCSP.y, oldCSP.chart.D0);
F = arl2EvalPsi(oldX, oldCSP.chart, data);
if abs(F - lm.J)>1.e-6*F,  F,  lm.J,  error('Bad lm'), end
%line 527 <arl2.nw>
% New CSP (Potapov)
allpass = arl2CSP2Allpass(oldCSP);
CSP     = arl2Allpass2CSP(allpass); 
% Force D0 = I in Potapov chart (harmless)
CSP.chart.D0 = eye(size(CSP.chart.D0));
chart   = CSP.chart;
p       = size(chart.D0,1);
%line 537 <arl2.nw>
% Add random unitary p-vector u and w=0 to chart
newChart.D0 = chart.D0;
newU = arl2CreateMatrix(p,1,'real'); 
newU = newU/norm(newU); 
newChart.u = [newU  chart.u ];
newChart.w = [ 0  chart.w  ];
%line 546 <arl2.nw>
newY = newU;
for k = 1:maxIP
 % Compatible Schur parameter
 rho = 0.8; rhomax = 0.999; tol = 1.e-6;
 while 1
   lambda = exp(2*i*(k-1)*pi/maxIP);
   Y = [lambda*rho*newY  CSP.y ];
   nX = arl2SP2X(Y, chart.D0);
   F = arl2EvalPsi(nX, newChart, data) ;
   if F <= (1+tol)*lm.J, break, end
   if rho>rhomax, error('Dichotomy failed'), end
   rho = (rhomax + tol + rho)/2;
 end
 newCSP.D = CSP.D;
 newCSP.y = Y;
 newCSP.chart=newChart;
 ipList{k} = arl2CSP2Allpass(newCSP);
end


