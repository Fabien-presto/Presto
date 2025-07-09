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
function plotb(varargin)
clf;
hold on;
k=2;
cmap=get(gcf,'DefaultAxesColorOrder');
[i,j]=size(cmap);
c=0;
hold on;
x=varargin{1};
while(k<=length(varargin))
    xy_aux=squeeze(varargin{k});
    y = 20*log10(abs(reshape(xy_aux,[length(xy_aux),1])));
    if ((k<length(varargin)) & isa(varargin{k+1},'char'))
        if (length(varargin{k+1})==1 & ~isletter(varargin{k+1}))
            plot(x,y,'color',cmap(c+1,:),'marker',varargin{k+1},'linestyle','none');
        else
            plot(x,y,varargin{k+1});
        end
        k=k+1;
    else
        plot(x,y,'color',cmap(c+1,:));
    end
    c=mod(c+1,i);
    k=k+1;
end
hold off; 

