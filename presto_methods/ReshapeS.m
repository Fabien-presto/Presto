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
function [S]=ReshapeS(S,check_flag)

if (nargin < 2)
	check_flag=1;
end
n=length(S.freq);
s=size(S.value);
assert(n==s(3),'S.value should have same length as S.freq');
if (check_flag ~=-1)
	assert(min(S.freq)*max(S.freq)<0,'S.freq should have sorted frequencies around zero - S is supposed to be low pass');
end

[S.freq,ind_s]=sort(S.freq);
S.value=S.value(:,:,ind_s);
S.freq=reshape(S.freq,1,n);

 

