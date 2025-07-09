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
function PlotS(S,graphic_type_flag,legend_flag)
% PlotS(S,graphic_type_flag,legend_flag)
% S: S informations to be ploted
% graphic_type_flag (optional): indicates if a Nyquist or a Bode plot is expected {'n','b'}
% legend_flag (optional): indicates if a legend should appear for each curve (data, compl., rat. app.) {0,1}
% Action
% Plots an S structure according to what it currently contains. For ex. if a completion is present
% it will be ploted in addition to the data. 

S=ReshapeS(S,-1);
DVC=GetDefaultValuesAndConstants(S);
[imax,jmax]=GetSDataSize(S);

% Standard value for plot
max_w_for_nyq=DVC.PS.max_w_for_nyq;
supp_ratio_for_bode=DVC.PS.supp_ratio_for_bode;

if (nargin<3)
    legend_flag=DVC.PS.legend_flag;
end
if (nargin<2)
   graphic_type_flag=DVC.PS.graphic_type_flag;
end

% Print name
if (isfield(S,'name'))
    set(gcf,'Name',S.name);
else
    set(gcf,'Name','');
end

clf;
subplot(imax,jmax,1);
if (graphic_type_flag=='b')
    PlotS_B(S,supp_ratio_for_bode,legend_flag);
else
    PlotS_N(S,max_w_for_nyq,legend_flag);
end

% Write subtitles
title_string=cell(imax,jmax);
for k=[1:imax]
    for j=[1:jmax]
        if (isfield(S,'comp'))
            if (isfield(S,'ratio'))
            st=sprintf(' r=%3.2f',S.ratio(k,j)*100);
            st=[st,'%'];
            title_string{k,j}=[title_string{k,j},st];
            end
        end
        if (isfield(S,'p') & isfield(S,'q'))
            Sid_proper_p=FreqResp(S.p,S.q,S.freq);
            vp=squeeze(Sid_proper_p(k,j,:));
            err=SimpleQuad(S.freq,abs(vp-squeeze(S.value(k,j,:))).^2);
            err=err/SimpleQuad(S.freq,abs(S.value(k,j,:)).^2);
            err=sqrt(err)*100;
            %err=norm(vp-squeeze(S.value(k,j,:)))/norm(squeeze(S.value(k,j,:)))*100;
            err_sup=norm(vp-squeeze(S.value(k,j,:)),inf)/norm(squeeze(S.value(k,j,:)),inf)*100;
            st=sprintf(' err=%3.2f',err);
            st=[st,'%'];
            st=[st,sprintf(' err.sup =%3.2f',err_sup)];
            st=[st,'%'];
            title_string{k,j}=[title_string{k,j},st];
        end    
        if (k~=j & imax==jmax)
            st=[];
            [diff,ind]=max(abs(S.value(1,2,:)-S.value(2,1,:)));
            max_12=(max(abs(S.value(1,2,:)))+max(abs(S.value(2,1,:))))/2;
            st=[st,sprintf(' diff-12-21 =%3.2f',100*diff/max_12)];
            st=[st,'%'];
            title_string{k,j}=[title_string{k,j},st];
        end
    end
end

for k=[1:imax]
    for j=[1:jmax]
        subplot(imax,jmax,jmax*(k-1)+j);
        hold on;
        title(title_string{k,j});
        hold off;
    end
end
 

