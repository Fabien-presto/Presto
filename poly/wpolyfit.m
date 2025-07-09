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
function [p,err]=wpolyfit(x,y,w,n)
m=length(x);
s=size(x);
if (s(1)==1)
    x=transpose(x);
end
s=size(y);
if (s(1)==1)
    y=transpose(y);
end
s=size(w);
if (s(1)==1)
    w=transpose(w);
end
y=w.*y;
A=w;
for j=2:n+1
    w=x.*w;
    A=[w,A];
end
AA=transpose(conj(A))*A;
yy=transpose(conj(A))*y;
xo=(AA)\(yy);
p=transpose(xo);
err=norm(A*xo-y); 

