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
%line 96 <utilities.nw>
function sys  = arl2CreateStableSystem(n,m,p,type)
if nargin<3, error('Syntax: sys  = arl2CreateStableSystem(n,m,p,type)'); end
if nargin < 4, type = 'complex'; end
   

if p<=m
 allpass  = arl2CreateAllpass(n,p,type);
 [A,B,C,D] = SSdata(allpass);
 B = arl2CreateMatrix(n,m,type);
else
 allpass  = arl2CreateAllpass(n,m,type);
 [A,B,C,D] = SSdata(allpass);
 C = arl2CreateMatrix(p,n,type);
end 
D = arl2CreateMatrix(p,m,type);
sys = SScreate(A,B,C,D);
