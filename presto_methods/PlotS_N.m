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
function PlotS_N(S,max_w,legend_flag)
[imax,jmax]=GetSDataSize(S);
legend_string=[];
for k=[1:imax]
    for j=[1:jmax]
        subplot(imax,jmax,jmax*(k-1)+j);
        v=squeeze(S.value(k,j,:));
        if (isfield(S,'p'))
            plot(real(v),imag(v),'.k');
            axis equal;
        else
            plot(real(v),imag(v),'k');
            axis equal;
        end
        if (k==1 & j==1 )
            legend_string=[legend_string;'data  '];
        end
        if (isfield(S,'comp'))
            if(k==1 & j==1)
               legend_string=[legend_string;'comp. ']; 
            end
            p=S.comp{k,j};
            w_c=[linspace((max(S.freq)),(max_w),10000),linspace((-max_w),(min(S.freq)),10000)];
            vv=polyval(p,1./w_c);
            hold on;
            if (isfield(S,'p') & isfield(S,'q'))
                plot(real(vv),imag(vv),'r','LineWidth',2);
            else
                plot(real(vv),imag(vv),'r');
            end
            hold off;
        end
    end
end
if (isfield(S,'p') & isfield(S,'q'))
    legend_string=[legend_string;'ratio.']; 
    theta=linspace(0,2*pi,1000);
    theta(1)=[];
    theta(end)=[];
    if (IsHighPass(S))
        w_eval=squeeze(S.freq);
    else
        w_eval=C2L(theta);
    end
    Sid_proper=FreqResp(S.p,S.q,w_eval);
    for k=[1:imax]
        for j=[1:jmax]
            v=squeeze(Sid_proper(k,j,:));            
            subplot(imax,jmax,jmax*(k-1)+j);
            hold on;
            plot(real(v),imag(v),'b');
            hold off;
        end
    end
end
subplot(imax,jmax,1);
if (legend_flag)
    legend(legend_string);
end
 

