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
function Matrix2Latex(A)
s=size(A);
fprintf('array(1..%d,1..%d,%s ',s(1),s(2),'[');
for i=1:s(1)
    fprintf(' %s ','[');
    for j=1:s(2)-1
        fprintf('%1.4f %s ',A(i,j),',');
    end
    fprintf('%1.4f %s ',A(i,j+1),']');
    if (i<s(1)) 
       fprintf(',\n');
   else
       fprintf('\n');
   end
end
fprintf(' %s \n ',']);');
 

