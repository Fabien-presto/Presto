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
function S=ApplyDiagonalRotation(S,a)

[imax,jmax]=GetSDataSize(S);

for k=1:imax
    for j=1:jmax
        if (k==j)
            S.value(k,k,:)=exp(i*a(k)).*S.value(k,k,:);
            if (isfield(S,'p'))
                S.p{k,k}=exp(i*a(k)).*S.p{k,k};
            end
            S.b(k,k)=angle(exp(i*S.b(k,k))*exp(i*a(k)));
        else
            if (k<=jmax)
                S.value(k,j,:)=exp(i*a(k)/2).*S.value(k,j,:);
                if (isfield(S,'p'))
                    S.p{k,j}=exp(i*a(k)/2).*S.p{k,j};
                end
                S.b(k,j)=angle(exp(i*S.b(k,j))*exp(i*a(k)/2));
            end
            if (j<=imax)
                S.value(k,j,:)=exp(i*a(j)/2).*S.value(k,j,:);
                if (isfield(S,'p'))
                    S.p{k,j}=exp(i*a(j)/2).*S.p{k,j};
                end
                S.b(k,j)=angle(exp(i*S.b(k,j))*exp(i*a(j)/2));
            end
        end
    end
end

    