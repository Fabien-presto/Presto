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
function [p,delay,errmin]=DetDelay(omega_in,y_in,delay_range,poly_order,zeros_at_inf,omega_lim,opt_flag,plot_flag)

if (nargin<8)
    plot_flag=0;
end
if (nargin<7)
    opt_flag=1;
end
if (nargin<6)
    zeros_at_inf=0;
end
l=length(omega_in);
assert(l==length(y_in),'omega and y should have same length');
omega=reshape(squeeze(omega_in),l,1);
y=reshape(squeeze(y_in),l,1);
[ind]=union(find(omega>abs(omega_lim)),find(omega<-abs(omega_lim)));
omega_inf=omega(ind);
y_inf=y(ind);
err=[];
if (length(omega_inf)<2*poly_order+1)
    fprintf('Warning: approximation degree is greater than 2x number of data points \n')
    assert(length(omega_inf)>poly_order,'Approximation degree is greater than number of data points');
end
%for r=delay_range
    %p=polyfit(1./omega_inf,exp(r*i.*omega_inf).*y_inf,poly_order);
 %   p=MyPolyfit(1./omega_inf,exp(r*i.*omega_inf).*y_inf,poly_order,zeros_at_inf);
  %  err=[err,norm(polyval(p,1./omega_inf)-((exp(r*i.*omega_inf)).*y_inf))];
  %end
[delay,p,err]=ComputeBestCompensation(omega_inf,y_inf,delay_range,poly_order,zeros_at_inf); 
if (plot_flag)
    plot(delay_range,err);
end
[errmin,indm]=min(err);
if ((indm==1 | indm==length(err)) & opt_flag)
    fprintf('Warning no local minimum found in delay det. - increases delay range \n');
end
step=sqrt((max(abs(omega_inf))-abs(omega_lim))/length(omega_inf));
errmin=errmin*step;
delay=delay_range(indm);
%p=MyPolyfit(1./omega_inf,exp(delay*i.*omega_inf).*y_inf,poly_order,zeros_at_inf);
 

