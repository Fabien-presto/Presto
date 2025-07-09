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
function WriteS(S,filename)
%  WriteS(S,filename)
% S: S-structure to write
% filename: destination filename
% Action
% Write the fields freq and value of S to an ascii file. Each line in
% the ascii file is structured as follows:
% S.freq(i) real(S.value(1,1,i)) imag(S(1,1,i)) real(S(1,2,i)) imag(S(1,2,i)) real(S(2,1,i)) imag(S(2,1,i) real(S(2,2i)) imag(S(2,2,i))

n=length(S.freq);
fid = fopen(filename,'w');
if (fid)
	for k=1:n
		fprintf(fid,'%1.9e ',S.freq(k));
		for m=1:2
			for l=1:2
				fprintf(fid,'%1.9e %1.9e ',real(S.value(m,l,k)),imag(S.value(m,l,k)));
			end
		end
		fprintf(fid,'\n');
	end
	fclose(fid);
else
	error('Cant write on file');
end
    



