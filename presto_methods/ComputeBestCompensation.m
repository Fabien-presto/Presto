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
function [d_opt,p,errs]=ComputeBestCompensation(w,y,d_range,poly_order,zeros_at_inf)
% Computes best delay compensation. Basically an optimized and optimized 
% least square solver for this problem

w=reshape(w,length(w),1);
y=reshape(y,length(y),1);

if (zeros_at_inf==0)
    P=GeneratePowerMatrix(1./w,poly_order);
else
    P=GeneratePowerMatrix(1./w,poly_order+zeros_at_inf);
    % Leave at least 2 freedom parameter
    i_max=max(poly_order+1,zeros_at_inf+3);
    P=P(:,[zeros_at_inf+1:i_max]);
    %P=P(:,[zeros_at_inf+1:poly_order+zeros_at_inf+1]);
end
tP=transpose(conj(P));
tPP=tP*P;
itPP=inv(tPP);
d_opt=-1;
e_opt=0;
errs=[];
for r=d_range
    yy=y.*exp(i*r*w);
    tPy=tP*(yy);
    v=itPP*tPy;
    err=norm(P*v-yy);
    errs=[errs,err];
    if (d_opt==-1 | err<e_opt)
        d_opt=r;
        e_opt=err;
        p=v;
    end
end  
p=[zeros([zeros_at_inf,1]);p];
p=transpose(p([length(p):-1:1]));

 

