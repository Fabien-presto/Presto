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
function [beta_opt,err_opt]=DetBeta(S,plot_flag)
n_beta_test=1000;
if (nargin<2)
    plot_flag=0;
end       
S11=squeeze(S.value(1,1,:));
S12=squeeze(S.value(1,2,:));
S21=squeeze(S.value(2,1,:));
S22=squeeze(S.value(2,2,:));
l=length(S22);
beta_range=linspace(0,pi,n_beta_test);
err=[];
for beta=beta_range
    err=[err,norm(S11.*conj(S12)+exp(i*2*beta)*S21.*conj(S22))/sqrt(l)];
end
if (plot_flag)
    plot(beta_range,err);
end
[err_opt,ind]=min(err);
beta_opt=beta_range(ind); 

