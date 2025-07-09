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
%line 115 <utilities.nw>
function allpass = arl2AllpassFactor(sys)
[n,m,p]   = SSdim(sys);
[A,B,C,D] = SSdata(sys);

if p<=m
 W = dlyap(A', C'*C);
 T = sqrtm(W);
 C = C/T;
 A = T * A / T;
 AC = [A;C];
 Z = AC*AC'; Z = eye(size(Z)) - Z;
 [U,S,V] = svd(Z);
 B = U(1:n,1:p);
 D = U(n+1:end,1:p);
else
 W = dlyap(A, B*B');
 T = sqrtm(W);
 B = T \ B;
 A = T \ A * T;
 AB = [A B];
 Z = AB'*AB; Z = eye(size(Z)) - Z;
 [U,S,V] = svd(Z);
 auxV=V';
 C = auxV(1:m,1:n); 
 D = auxV(1:m,n+1:end); 
end

allpass = SScreate(A,B,C,D);
arl2Assert(arl2IsAllpass(allpass));
