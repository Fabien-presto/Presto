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
function [H]=GenerateMoebMatrix(theta1,theta2,n_cf,n)

n_sample=2000;
theta=linspace(theta1,theta2,n_sample);
z=exp(i*theta);
H=[];
for m=[0:n-1]
	if (m==0)
		v=[];
		for k=[1:n_cf]
			v=[v;(exp(-i*k*theta2)-exp(-i*k*theta1))/(-i*k*2*pi)];
		end
		H=[v];
	else
		H=[H,transpose(FourierCoeffs(theta,(i*moeb(z,1)).^m,[1:n_cf]))];
	end
end	



 

