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
function [delay,svdmin]=ComputeBestCausalCompensation(w,y,d_range,poly_order,deg_rational,n_cf,w_lim,modulus_factor)
% Computes best delay compensation by looking for best causal completion 
% delay is the optimal delay found and ecmin is a measure of the non causal part

number_of_spline_points=300;
ws=linspace(min(w),max(w),number_of_spline_points);
ys=spline(w,y,ws);
%ws=w;
%ys=y;
k=0;
e=zeros(size(d_range));
for r=d_range
    k=k+1;
    yc=ys.*exp(i*r*ws);
    [p,q,e(k)]=l2ratappeig(ws,yc,deg_rational,deg_rational);
end
clf;
plot(d_range,e);
pause;
ipos=find(d_range>0);
[svdmin,i_min]=min(e(ipos))
d_range_p=d_range(ipos);
delay=d_range_p(i_min);
