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
function [SS,sys]=RatApp(S,n,divide_by_z_minus_one,solver_flag); 
% [Sres]=RatApp(S,n,divide_by_z_minus_one,solver_flag);
% S: S measurements
% n: target MacMillan degree
% divide_by_z_minus_one: fixes the type of isometry used to pass to the circle - 
%                        poss. values: {0,1} - 0 is L^{\infty} isom. 1 is is L^2 isom.                            
% solver_flag: specifies which rational approximation engine to use
% between Rarl2 and hyperion. poss.value: {'RARL2','HARL2'}
% Action
% Launches a rational approximation engine for S
% Important: Completion and Fourier coefficients should be computed before 
%            see: CompensateDelayAndFreqShift, CompFourierCoeffs
% The result is S enriched the rational approximation (S.p and S.q)


assert(isfield(S,'fourier'),'Compute fourier coeffs first');
DVC=GetDefaultValuesAndConstants(S);
sign_in_zero_trans_12=DVC.CDAFS.sign_in_zero_trans_12;
adjust=DVC.RA.shift_adjustment & IsFilter(S);
adjust_post=DVC.RA.shift_adjustment_post & IsFilter(S);
stop_tresh=DVC.RA.func_stop_crit_optim;
stop_tresh_post=DVC.RA.func_stop_crit_post_optim;
[imax,jmax]=GetSDataSize(S);
if (nargin<4)
    solver_flag=DVC.RA.solver_flag;
end
if (nargin<3)
    divide_by_z_minus_one=DVC.RA.divide_by_z_minus_one;
end
SS=S;
% Divide by z-1 if necessary
if (divide_by_z_minus_one)    
            [S.fourier,vals_in_one]=DivideByZminusOne(S.fourier);
end

post_optimization=DVC.RA.post_optimization;

f_proper=cell(imax,jmax);
% Take proper part 
for k=1:imax
    for j=1:jmax
        f_proper{k,j}=(S.fourier{k,j}(2:end));
    end
end
% Launch solver
switch (solver_flag)
    case 'HARL2'
        [p,q]=Run_Harl2(f_proper,n);
    case 'MRARL2'
        tic;
        [p,q,sys,beta]=Run_MRarl2(f_proper,n,adjust,stop_tresh);
        fprintf('Time to compute core-optim  in rat. app. %f\n',toc);
        if (adjust)
            if (divide_by_z_minus_one)
                [SS,sys,p]=ApplyAdditionalBeta(SS,beta,sign_in_zero_trans_12,sys,p,vals_in_one);
            else
                [SS,sys,p]=ApplyAdditionalBeta(SS,beta,sign_in_zero_trans_12,sys,p,S.fourier);
            end
        end
    case 'AAK'
        tic;
        [p,q,sys,beta]=Run_AAK(f_proper,n,adjust,stop_tresh);
        fprintf('Time to compute core-optim by AAK only %f\n',toc);
    otherwise
        tic;
        [p,q,sys,beta]=Run_Rarl2(f_proper,n,adjust,stop_tresh);
        fprintf('Time to compute core-optim  in rat. app. %f\n',toc);
        if (adjust)
            if (divide_by_z_minus_one)
                [SS,sys,p]=ApplyAdditionalBeta(SS,beta,sign_in_zero_trans_12,sys,p,vals_in_one);
            else
                [SS,sys,p]=ApplyAdditionalBeta(SS,beta,sign_in_zero_trans_12,sys,p,S.fourier);
            end
        end
end
% Adding proper part
for k=1:imax
    for j=1:jmax
        if (divide_by_z_minus_one)
            if (k~=j)
                vals_in_one(k,j)=0;
            end
            p{k,j}=polyadd(polymult(p{k,j},[1 -1]),vals_in_one(k,j)*q);
            if (post_optimization)
                sys.d(k,j)=S.fourier{k,j}(1);
            end
        else
            p{k,j}=polyadd(p{k,j},S.fourier{k,j}(1)*q);
            if (post_optimization)
                sys.d(k,j)=S.fourier{k,j}(1);
            end
        end
    end
end
if (divide_by_z_minus_one)
    sys.d=sys.d+sys.c*sys.b;
    sys.c=sys.c*sys.a-sys.c;
end

% Launch poste l2 optimization
if (post_optimization)
    thetas=L2C(S.freq);
    switch(solver_flag)
        case 'MRARL2'
            tic;
            [p,q,beta]=Run_MRarl2_l2(thetas,SS.value,sys,n,adjust_post,stop_tresh_post);
            fprintf('Time to compute post-optim  in rat. app. %f\n',toc);
            if (adjust_post)
                [SS,sys,p]=ApplyAdditionalBeta(SS,beta,sign_in_zero_trans_12,sys,p);
            end
        otherwise
            tic;
            [p,q,beta,sys]=Run_Rarl2_l2(thetas,SS.value,sys,n,adjust_post,stop_tresh_post);
            fprintf('Time to compute post-optim  in rat. app. %f\n',toc);
            if (adjust_post)
                [SS,sys,p]=ApplyAdditionalBeta(SS,beta,sign_in_zero_trans_12,sys,p);
            end
    end
end
% Passing to plane
[p,q]=MoebSubsInTrans(p,q);

% Clean simplification if the diagonal are supposed to be theoreticaly strictly proper
if (divide_by_z_minus_one & ~post_optimization)
    for k=1:imax
        for j=1:jmax
            if (k~=j & length(p{k,j})==n+1)
                p{k,j}(1)=[];
            end
        end
    end
end
SS.p=p;
SS.q=q;
 

