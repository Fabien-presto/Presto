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
%line 111 <allpass.nw>
function [allpass, Slist, Ulist, Vlist] = arl2CSP2Allpass(CSP)

% Need gradients ?
switch nargout
 case 1, needGrad = 0;
 case 4, needGrad = 1;
 otherwise, error('Bad output args count')
end

% Dimension and degree
[p, N] = size(CSP.y);

% First 0-order allpass system
DCBA = CSP.chart.D0;
%line 127 <allpass.nw>
% Loop on degree
for n = 1:N
%line 131 <allpass.nw>
 % Build U,V matrices
 [U, V] = arl2UV(CSP, n);
 Id=eye(n-1);
 Z = zeros(p+1,n-1);
 UU = [ U' Z ; Z' Id ];
 VV = [ V  Z ; Z' Id ]; 
%line 139 <allpass.nw>
 % Apply (FUV) transform
 Z = zeros(1,p+n-1);
 S = [ 1 , Z ;  Z', DCBA ];
 DCBA = VV * S * UU;
 arl2Assert(arl2IsUnitary(DCBA));
%line 146 <allpass.nw>
 if needGrad
  Ulist{n} = UU;
  Vlist{n} = VV;
  Slist{n} = S;
 end
%line 153 <allpass.nw>
% End loop on degree
end
%line 157 <allpass.nw>
% Buils allpass
D = DCBA(1:p,  1:p);
C = DCBA(1:p,  p+1:end);
B = DCBA(p+1:end, 1:p);
A = DCBA(p+1:end, p+1:end);
allpass = SScreate(A,B,C,D);
%line 166 <allpass.nw>
function [U, V] = arl2UV(CSP, n)

u = CSP.chart.u(:,n);
w = CSP.chart.w(n);
v = CSP.y(:,n);

Iw2=sqrt(1-w'*w);
Iv2=sqrt(1-v'*v);
Iwv2=sqrt(1-w'*w*v'*v);
uu=u*u';
Id=eye(size(uu));

U = [ (Iw2/Iwv2)*u  ,  Id-(1 + w*Iv2/Iwv2)*uu
      w'*Iv2/Iwv2   ,  (Iw2/Iwv2)*u'                  ];

if norm(v) > 0
V = [ (Iw2/Iwv2)*v  ,  Id-(1 - Iv2/Iwv2)*v*v'/(v'*v)
      Iv2/Iwv2      ,  -(Iw2/Iwv2)*v'                 ];
else
V = [ (Iw2/Iwv2)*v  ,  Id
      Iv2/Iwv2      ,  -(Iw2/Iwv2)*v'                 ];
end

% Check
arl2Assert(arl2IsUnitary(U));
if ~arl2IsUnitary(V)
 disp(V)
 disp(norm (eye(size(V)) - V'*V))
 dlbquit
end

