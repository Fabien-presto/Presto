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
%line 301 <arl2.nw>
function [C , Ceq, Cgrad, Ceqgrad] = arl2SPModulus(X, chart, varargin)
[p,n] = size(chart.u);
Z = arl2X2SP(X,n,p);
if p>1, C = sum(abs(Z).^2); else, C = abs(Z).^2 ;end

Rmax = 0.95;
C = C  - Rmax;

Cgrad = [];
for constr = 1:n
 ZZ = zeros(size(Z)); 
 ZZ(:,constr) = Z(:,constr);
 ZZ = ZZ(:);
 Cgrad = [Cgrad , [real(ZZ);imag(ZZ)]];
end
d= length(X)-2*n*p;
Cgrad = [Cgrad; zeros(d, size(Cgrad,2))];
Ceqgrad = [];
Ceq = [];
