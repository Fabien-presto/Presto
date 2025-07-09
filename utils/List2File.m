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
function List2File(filename,LL)

fp=fopen(filename,'w');
for i=1:length(LL)
	n=length(LL{i});
	A=LL{i};
	fprintf(fp,'Sol(%d): \n',i);
	for k=1:n
		for l=1:n
			if (A(k,l)<0) 
			  fprintf(fp,' %1.4f',real(A(k,l)));
			else 
                if (A(k,l)==0)
                    fprintf(fp,'   0   ');
                else
		  fprintf(fp,'  %1.4f',abs(real(A(k,l))));
                end
			end
		end
		fprintf(fp,'\n');
	end;
	fprintf(fp,'\n\n');
end;
fclose(fp);
