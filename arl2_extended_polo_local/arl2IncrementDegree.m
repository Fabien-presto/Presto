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
%line 352 <arl2.nw>
function [newX, newChart]  = arl2IncrementDegree(oldX, oldChart)

% Build current allpass system
allpass = arl2X2Allpass(oldX, oldChart);

% Recode it
newCSP   = arl2Allpass2CSP(allpass);
newChart = newCSP.chart;

% Add u and w
[p, n] = size(newChart.u);
w = 0; u = arl2CreateMatrix(p,1); u = u/norm(u);
newChart.u = [u newChart.u];
newChart.w = [w newChart.w];
newCSP.chart = newChart;

% Add y
y = newChart.D0*u;
newCSP.y = [y newCSP.y];
newX     = arl2SP2X(newCSP.y, newCSP.D);

% Recode !
%[newX, newChart]  = arl2ChangeChart(newX, newChart);
