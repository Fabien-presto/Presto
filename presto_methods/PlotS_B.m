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
function PlotS_B(S,s_ratio,legend_flag)
[imax,jmax]=GetSDataSize(S);
legend_string=[];
h=[];
for k=[1:imax]
    for j=[1:jmax]
        subplot(imax,jmax,jmax*(k-1)+j);
        v=20*log10(abs(squeeze(S.value(k,j,:))));
        if (isfield(S,'p'))
            h_aux=plot(S.freq,v,'.k');
        else
            h_aux=plot(S.freq,v,'k');
        end
        if (k==1 & j==1 )
            legend_string=[legend_string;'data  '];
            h=[h;h_aux];
        end
        if (isfield(S,'comp'))
            p=S.comp{k,j};
            w1=linspace(max(S.freq),max(S.freq)*s_ratio,500);
            w2=linspace(min(S.freq)*s_ratio,min(S.freq),500);
            vv1=20*log10(abs(polyval(p,1./w1)));
            vv2=20*log10(abs(polyval(p,1./w2)));
            hold on;
            if (isfield(S,'p') & isfield(S,'q'))
                h_aux=plot(w1,vv1,'r','LineWidth',2);
                plot(w2,vv2,'r','LineWidth',2);
            else
                h_aux=plot(w1,vv1,'r');
                plot(w2,vv2,'r');
            end
            if(k==1 & j==1)
               legend_string=[legend_string;'comp. '];
               h=[h;h_aux];
            end
            hold off;
        end
    end
end
if (isfield(S,'p') & isfield(S,'q'))
    if (IsHighPass(S))
        w_eval=linspace(min(S.freq)/s_ratio,max(S.freq)*s_ratio,1000);
    else
        w_eval=linspace(min(S.freq)*s_ratio,max(S.freq)*s_ratio,1000);
    end
    Sid_proper=FreqResp(S.p,S.q,w_eval);
    for k=[1:imax]
        for j=[1:jmax]
            v=20*log10(abs(squeeze(Sid_proper(k,j,:))));            
            subplot(imax,jmax,jmax*(k-1)+j);
            hold on;
            h_aux=plot(w_eval,v,'b');
            hold off;
            if(k==1 & j==1)
               legend_string=[legend_string;'ratio.'];
               h=[h;h_aux];
            end
        end
    end
end
subplot(imax,jmax,1);
if (legend_flag)
    legend(h,legend_string);
end
 

