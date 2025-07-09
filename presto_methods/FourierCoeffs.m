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
function [f]=FourierCoeffs(theta,y,k_range)
theta=reshape(theta,1,length(theta));
y=reshape(y,1,length(y));
l=length(theta);
eps0=(theta(2:l)-theta(1:l-1))/2;
eps1= [eps0,0];
eps2= [0,eps0];
eps=(eps1+eps2);
yw =eps.*y;
f=[];
v0=exp(-k_range(1)*i*theta);
vd=exp(-i*theta);
for k=k_range
    %c=sum(exp(-k*i*theta).*yw)/2/pi;
    c=sum(v0.*yw)/2/pi;
    f=[f,c];
    v0=v0.*vd;
end
 

