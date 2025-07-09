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
%line 33 <userF.nw>
function data = arl2Preprocess(userData)
data = userData;
switch data.type
case 'sys'
 [a,b,c,d] = SSdata(userData.value);
 data.D = d;
 data.F0 = SSnH2(data.value);
case 'sample'
 data.D  = mean(data.value,3);
 data.F0 = real(sum(sum(mean(data.value .* conj(data.value) , 3))));
case 'coef'
 data.D  = data.valinf;
 x = data.D(:); y = data.value(:);
 data.F0 = real(x'*x + y'*y);
end
