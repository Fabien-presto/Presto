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
function [p]=MyPolyfit(x,y,n,n_prescibed_zero_in_zero)
if (n_prescibed_zero_in_zero==0)
    p=polyfit(x,y,n);
else
    assert(n>=n_prescibed_zero_in_zero,'Prescribed zero at 0 should be less or equal total degree');
    % Build matrix
    nx=length(x);
    y=reshape(y,nx,1);
    A=zeros(nx,n+1-n_prescibed_zero_in_zero);
    for k=1:nx
        A(k,n-n_prescibed_zero_in_zero+1)=x(k)^(n_prescibed_zero_in_zero);
        for l=(n-n_prescibed_zero_in_zero):-1:1
            A(k,l)=A(k,l+1)*x(k);
        end
    end
    tA=transpose(conj(A));
    p=(tA*A)\(tA*y);
    p=[transpose(p),zeros(1,n_prescibed_zero_in_zero)];
end 

