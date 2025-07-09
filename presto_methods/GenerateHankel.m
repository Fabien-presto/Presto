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
function [H]=GenerateHankel(theta1,theta2,n,m)
% Generates the Hankel matrix associated with the linear operator
% P_{\bar H^2}(g.Chi_J) where J is the interval the subarc of the unit
% circle determined by theta1 and theta2 (trigo sens). n,m is the size
% of the resulting H     
  
H=zeros(n,m);
for k=1:n
    H(k,1)=(exp(i*k*theta2)-exp(i*k*theta1))/(k*2*pi*i);
    for j=1:min(m-1,k-1)
        H(k-j,j+1)=H(k,1);
    end
end
for l=2:m
    k=n+l-1;
    H(n,l)=(exp(i*k*theta2)-exp(i*k*theta1))/(k*2*pi*i);
    for j=l+1:m
        H(n+l-j,j)=H(n,l);
    end
end
 

