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
function C=MatrixProd3D(A,B)
% Computes the matricial product of the 3d matrices A an B along their third dimension  
s1=size(A);
s2=size(B);

assert(length(s1)==3 && length(s2)==3,'Matrices should have 3 diemsnions');
assert(s1(3)==s2(3),'Matrices should have same depth');
assert(s1(2)==s2(1),'Matrices dimension mismatch for mulitplication');
if (s1(1)>s2(2))
    C=permute(MatrixProd3D(permute(B,[2,1,3]),permute(A,[2,1,3])),[2,1,3]);
else
    C=zeros(s1(1),s2(2),s2(3));
    for k=1:s1(1)
        C(k,:,:)=sum(repmat(permute(A(k,:,:),[2,1,3]),[1,s2(2)]).*B,1);
    end
end

