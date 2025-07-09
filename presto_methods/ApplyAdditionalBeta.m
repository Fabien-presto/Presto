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
function [S,sys,p]=ApplyAdditionalBeta(S,beta,signe_rule,sys,p,val_anti_diag)

[x,k_zero]=min(abs(S.freq));
% Fix the sign of the imaginary part of S12 in zero 
if (IsFilter(S) & imag(exp(i*beta)*(S.value(1,2,k_zero)))*signe_rule<0 )
    beta=beta+pi;
    if (nargin >3) 
        sys.c(1,:)=-1*sys.c(1,:);
        sys.b(:,1)=-1*sys.b(:,1);
        p{1,2}=-p{1,2};
        p{2,1}=-p{2,1};
    end
end
m=exp(i*beta);
[imax,jmax]=GetSDataSize(S);
for k=1:imax
    for j=1:jmax
        if (k~=j)
            S.value(k,j,:)=m*S.value(k,j,:);
            if (isfield(S,'comp'))
                S.comp{k,j}=m*S.comp{k,j};
            end
        end
    end
end
if (nargin>6) 
    for k=1:imax
        for j=1:jmax
            if (k~=j)
                val_anti_diag(k,j)=m*val_anti_diag(k,j);
            end
        end
    end
end
S.b(1,2)=S.b(1,2)+beta;
S.b(2,1)=S.b(2,1)+beta;