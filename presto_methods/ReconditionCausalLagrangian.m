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
function [L,Pc]=ReconditionCausalLagrangian(L,x1,x2)
% Reconditonne le lagrangien en changeant de base pour la completion
% La nouvelle base est est constituee de poly. de Chebychev sur le segment [x1,x2]
% Pc est la matrice de changement de base

s=size(L.D);
n=s(2)-1;
Paux=transpose(ChebychevBasis(x1,x2,n));
Pc=Paux(end:-1:1,:);

L.D=L.D*Pc;
L.B=L.B*Pc;

if (isfield(L,'BB'))
    L.BB=L.B'*L.B;
end
