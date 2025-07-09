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
%line 243 <utilities.nw>
function SSplotMat(sample, varargin)
if ~isempty(varargin), color=varargin{1}; else, color = 'b'; end
global SSplotFigure
if ~isempty(SSplotFigure), figure(SSplotFigure); 
else, return, end
[p,m,L] = size(sample);
index = 0;
for lig=1:p
for col=1:m
 index = index + 1;
 subplot(p,m,index);
 resp=squeeze(sample(lig,col,:));
 plot(real(resp), imag(resp),color);
 axis equal
 hold on;
end
end
drawnow
SSplotFigure = gcf;
