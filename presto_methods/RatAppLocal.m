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
function [Sr,sys,flag]=RatAppLocal(S,sys)
% [Sr,sys,flag]=RatAppLocal(S,sys)
% Computes a rational model of the data present in S by using the bootstrap
% structure sys - the optimisation is done ine l2 terms directly on the data 

DVC=GetDefaultValuesAndConstants(S);
sign_in_zero_trans_12=DVC.CDAFS.sign_in_zero_trans_12;
adjust=DVC.RA.shift_adjustment & IsFilter(S);
adjust_post=DVC.RA.shift_adjustment_post & IsFilter(S);
stop_tresh=DVC.RA.func_stop_crit_optim;
stop_tresh_post=DVC.RA.func_stop_crit_post_optim;
solver_flag=DVC.RA.solver_flag;
n=length(sys.a);
flag='true';

thetas=L2C(S.freq);
try
    switch(solver_flag)
        case 'MRARL2'
            tic;
            [p,q,beta]=Run_MRarl2_l2(thetas,S.value,sys,n,adjust_post,stop_tresh_post);
            fprintf('Time to compute post-optim  in rat. app. %f\n',toc);
            if (adjust_post)
                [S,sys,p]=ApplyAdditionalBeta(S,beta,sign_in_zero_trans_12,sys,p);
            end
        otherwise
            tic;
            [p,q,beta,sys]=Run_Rarl2_l2(thetas,S.value,sys,n,adjust_post,stop_tresh_post);
            fprintf('Time to compute post-optim  in rat. app. %f\n',toc);
            if (adjust_post)
                [S,sys,p]=ApplyAdditionalBeta(S,beta,sign_in_zero_trans_12,sys,p);
            end
    end
catch err
    err
    flag='false';
    warning('Local Optimsation failed');
    Sr=S;
end
% Passing to plane
[p,q]=MoebSubsInTrans(p,q);
Sr=S;
Sr.p=p;
Sr.q=q;
Sr=AdjustValueAtInfinityFromRattApp(Sr);