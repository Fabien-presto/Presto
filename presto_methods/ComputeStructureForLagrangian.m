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
function [L]=ComputeStructureForLagrangian(x,y,xc,yc,fs,causalc,n,cflag,iso_flag)
% Builds up the structure used by EvaluateLagrangian for evaluation
% purposes
% This structure contains all the necessary parmeters to compute the
% value at a dual point lambda.

x=reshape(x,length(x),1);
y=reshape(y,length(y),1);
xc=reshape(xc,length(xc),1);
yc=reshape(yc,length(yc),1);
L.P=GeneratePowerMatrix(x,n);
L.Pc=GeneratePowerMatrix(xc,n);
L.y=y;
L.yc=yc;
L.f=reshape(fs.f,length(fs.f),1);
% Generate two different matrices for estimation of the non-causal part
% depending of the kind of isometry used.
if (iso_flag)
	L.H=GenerateDivMoebMatrix(fs.theta1,fs.theta2,length(L.f),n+1);
else	
	L.H=GenerateMoebMatrix(fs.theta1,fs.theta2,length(L.f),n+1);
end
L.causalc=causalc;
L.norme_f2=fs.norm_f2;

% Flag concerning the type of contrained problem to be solved
L.cflag=cflag;
 

