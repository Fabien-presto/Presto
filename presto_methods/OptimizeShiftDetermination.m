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
function [beta]=OptimizeShiftDetermination(fouriers,n)

m=100;
beta_range=linspace(-pi/50,pi/50,m);

crits=[];
mfouriers=fouriers;
for r=beta_range
    edelta=exp(i*r);
    for k=1:2*50+1
        fouriers{1,2}(k)=mfouriers{1,2}(k)*edelta;
        fouriers{2,1}(k)=mfouriers{2,1}(k)*edelta;
    end
    for k=1:2*50+1
        for l=1:2
            for m=1:2
                f(l,m)=fouriers{l,m}(k+1);
            end
        end
        blocs{k}=f;
    end
    H=CreateBlocHankel(blocs);
    sl=svds(H,n+1)
    crits=[crits,sl(n+1)];
end
[critmin,imin]=min(crits);
beta=beta_range(imin);
plot(beta_range,crits);
        