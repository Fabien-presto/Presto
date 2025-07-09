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
%line 39 <utilities.nw>
function matPrint(A)
[lig,col]=size(A);
if ~isempty(A)
 if isreal(A), type = 'real'; else, type = 'complex'; end
 for i=1:lig
  for j=1:col
   x = A(i,j);
   printItem(type, x)
  end
  fprintf('\n');
 end
end
fprintf('\n');
%line 54 <utilities.nw>
function printItem(type, x)
switch type
case 'complex'
 if imag(x)<0 , sig = ' -'; else, sig = ' +'; end
 fprintf(' %6.3f%2s%6.3fi ',real(x),sig,abs(imag(x)));
case 'real'
 fprintf(' %6.3f',x);
otherwise
 mesg = sprintf('%s: unknown type',type);
 error(mesg);
end
