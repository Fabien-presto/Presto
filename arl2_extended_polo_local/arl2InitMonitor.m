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
%line 402 <arl2.nw>
function arl2InitMonitor(f, n, xmax , nbDec)
global monitorData

monitorData.status = 'on';
H = 50*(n+1);
pos = [ 400, 10,   450,  H];

if ~isfield(monitorData,'figure') | isempty(monitorData.figure)
 monitorData.figure = figure('Position',pos);
else
 figure(monitorData.figure);
 set(gcf,'Position',pos);
end
set(gcf,'Name','Function & Schur parameters (module)');
clf
monitorData.fmax    = f;
monitorData.xmax   = xmax;
monitorData.shift  = 0;
if nargin < 4, monitorData.nbDec = 6; 
else, monitorData.nbDec = nbDec; end
