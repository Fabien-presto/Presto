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
%line 37 <allpass.nw>
function CSP = arl2Allpass2CSP(sys)

if ~arl2IsAllpass(sys), disp(sys), error('Not allpass system'), end
%line 42 <allpass.nw>
% Extract system matrices
[A,B,C,D] = SSdata(sys);
n = size(A,1);
p = size(D,1);

% Unitary conversion to lower triangular form
if n > 0
 [U,T] = schur(A','complex') ;
 B = U'*B;
 C = C*U;
 A = T';
end
%line 56 <allpass.nw>
% Set Schur parameters (adapted chart y==0)
CSP.y = zeros(p,n);

% Init chart
chart.u = [];
chart.w = [];
%line 64 <allpass.nw>
% Loop on degree
while n > 0
%line 68 <allpass.nw>
 % Split system
 index = 2:n;
 a = A(1,1);
 beta = A(index,1);
 b = B(1,:);
 c = C(:,1);
 A = A(index,index);
 B = B(index,:);
 C = C(:,index);
%line 79 <allpass.nw>
 % Build current chart parameters
 w = a;
 arl2Assert(arl2IsSmall(norm(b*b'+w*w'-1)));
 u = b'/norm(b);
%line 85 <allpass.nw>
 % Store chart parameters
 chart.u = [u chart.u];
 chart.w = [w chart.w];
%line 90 <allpass.nw>
 % Updata D, B
 x = (1+w')/(1-w*w');
 D = D + x*c*b;
 B = B + x*beta*b;
 arl2Assert(arl2IsUnitary([A,B;C,D]));
%line 97 <allpass.nw>
 % End loop on degree
 n = n - 1;
end
%line 102 <allpass.nw>
% Return CSP
chart.D0  = D;
% matPrint(D)
CSP.chart = chart;
CSP.D = eye(size(D));
