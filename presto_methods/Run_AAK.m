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
function [p,q,sys,beta]=Run_AAK(fourier_,n,adjust,stop_treshold)

% Put fourier coeffs in an array
max_size=0;
[imax,jmax]=size(fourier_);
for k=1:imax
    for j=1:jmax
        max_size=max(max_size,length(fourier_{k,j}));
    end
end
fourier=zeros(imax,jmax,max_size);
for k=1:imax
    for j=1:jmax
        l=length(fourier_{k,j});
        fourier(k,j,1:l)=fourier_{k,j};
    end
end

%sys0 = arl2CreateStableSystem(n,2,2);
sys=fc2ir(fourier,n);
beta=0;

% Build transfer from state space description

[num,den]=MyTf(sys.a,sys.b,sys.c,sys.d);
for k=1:imax
    for j=1:jmax
        if (length(num{k,j})==n+1)
            num{k,j}(1)=[];
        end
    end
end
p=num;
q=den{1,1}; 

