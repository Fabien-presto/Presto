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
function [P]=ChebychevBasis(x1,x2,n)
% Returns a Chebychev basis of chebychev polynmials of max order n on the
% segment [x1,x2]
% Stores polynomials lines by lines

P=zeros(n+1,n+1);

a=(x2-x1)/2;
b=(x2+x1)/2;

P(1,n+1)=1;
for k=1:n
    v=(pi/2+[0:1:k-1]*pi)/k;
    v=cos(v)*a+b;
    aux=poly(v);
    aux=aux/polyval(aux,a+b);
    P(k+1,[n+1-k:n+1])=aux;
end
    