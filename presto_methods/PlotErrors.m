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
function PlotErrors(S)

w=squeeze(S.freq);
assert(isfield(S,'p')& isfield(S,'q'),'No rational approximation attached to the S-structure');
sd=FreqResp(S.p,S.q,w);
subplot(2,2,1);
for k=1:2
    for j=1:2
        subplot(2,2,2*(k-1)+j);
        plot(S.freq,abs(squeeze(sd(k,j,:))-squeeze(S.value(k,j,:))));
        %plot(S.freq,abs(sd(k,j)));
        %plot(S.freq,abs(squeeze(S.value(k,j,:))));
    end
end
