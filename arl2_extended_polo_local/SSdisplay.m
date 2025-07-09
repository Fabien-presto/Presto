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
%line 265 <utilities.nw>
function SSdisplay(sys)

[A,B,C,D]=SSdata(sys);

if isreal([A B; C D])
   type = 'real'; 
   whiteCol = '      ';
else
   type = 'complex'; 
   whiteCol = '               ';
end

[ligA,colA]=size(A);
[ligB,colB]=size(B);
[ligC,colC]=size(C);
[ligD,colD]=size(D);

%line 284 <utilities.nw>
if ~isempty(A)
 fprintf('a ='); 
 for j=1:colA, fprintf(whiteCol); end, if ~isempty(B), fprintf('b =\n'); end
 for i=1:ligA
 for j=1:colA
  x = A(i,j);
  printItem(type, x)
 end
 fprintf('   ');
 for j=1:colB
  x = B(i,j);
  printItem(type, x)
 end
 fprintf('\n');
 end
 fprintf('\n');
end
%line 303 <utilities.nw>
if ~isempty(C), fprintf('c ='); for j=1:colC, fprintf(whiteCol); end, end
if ~isempty(D), fprintf('d =\n'); end
for i=1:ligC
 for j=1:colC
  x = C(i,j);
  printItem(type, x)
 end
 if ~isempty(C), fprintf('   '); end
 for j=1:colD
  x = D(i,j);
  printItem(type, x)
 end
 fprintf('\n');
end
fprintf('\n');

fprintf('\n');
%line 322 <utilities.nw>
function printItem(type, x)
switch type
case 'complex'
 if imag(x)<0 , sig = ' -'; else, sig = ' +'; end
 fprintf(' %5.2f%2s%5.2fi ',real(x),sig,abs(imag(x)));
case 'real'
 fprintf(' %5.2f',x);
otherwise
 mesg = sprintf('%s: unknown type',type);
 error(mesg);
end
