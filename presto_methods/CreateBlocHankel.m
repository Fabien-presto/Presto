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
function [H]=CreateBlocHankel(blocs)

n=length(blocs);
assert( n>0 & mod(n,1)==0,'No blocs or no odd number of blocs');
[l,m]=size(blocs{1});
nn=(n+1)/2;
H=zeros(l*nn,m*nn);
for k=1:nn
    for j=1:k
        assert(prod([l,m]==size(blocs{k})),'Blocs should be of same size'); 
        H((k-j)*l+1:(k-j+1)*l,(j-1)*m+1:(j)*m)=blocs{k};        
    end
end
offset=l*nn;
for k=nn+1:n
    kk=k-nn;
    for j=1:nn-kk
        assert(prod([l,m]==size(blocs{k})),'Blocs should be of same size'); 
        H(offset+(kk-j-1)*l+1:offset+(kk-j)*l,j*m+1:(j+1)*m)=blocs{k};        
    end
end