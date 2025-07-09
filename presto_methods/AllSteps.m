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
function [S,sys]=AllSteps(S,n,zeros_at_inf,solver_flag,plot_flag)
% [Sres,sys]=AllSteps(S,n,zeros_at_inf,solver_flag,plot_flag)
% S: S measurements
% n: target MacMillan degree
% zeros_at_inf (optional): annulation order at infty for S12 and S21
% solver_flag (optional) : rational approx. engine used - {'RARL2','HARL2'}
% plot_flag (optional) : {0,1}
% Action
% Proceeds to the following steps:
%          - CompensateDelayAndFreqShift (and build completion)
%          - CompFourierCoeffs
%          - RatApp (launches rational approximation engine)
% The result is S enriched with the completion, the compensations, the rational approximation
% The output sys is a bootstrap strucutre than can be used in further
% identifications as starting point - typically in the function RatAppLocal
 

assert(nargin>1,'No McMillan degree given');
S=ReshapeS(S);
% Get Default values and Constants
DVC=GetDefaultValuesAndConstants(S);

if (nargin<5)
    plot_flag=DVC.AS.plot_flag;
end
if (nargin<4 | isempty(solver_flag))
    solver_flag=DVC.AS.solver_flag;
end
if (nargin<3 | isempty(zeros_at_inf))
    zeros_at_inf=DVC.AS.zeros_at_inf;
end
if (nargin==3 & ischar(zeros_at_inf))
    solver_flag=DVC.AS.solver_flag;
    zeros_at_inf=DVC.AS.zeros_at_inf;
end
if (plot_flag)
    clf;
    PlotS(S);
    RaiseAndDraw;
end
fprintf('Compensate delay and freq. shift - compute completion \n');
tic;
S=CompensateDelayAndFreqShift(S,zeros_at_inf);
fprintf('Time to compensate delays and shift:%f\n',toc);
if (plot_flag)
    PlotS(S);
    RaiseAndDraw;
    %pause;
end
fprintf('Compute  fourier coeff. - take stable part\n');
tic;
S=CompFourierCoeffs(S);
fprintf('Time to compute fourier coeffs:%f\n',toc);
if (plot_flag)
    PlotS(S);
    RaiseAndDraw;
end
fprintf('Compute  rational approximation \n');
t_rat_app=tic();
[S,sys]=RatApp(S,n,DVC.RA.divide_by_z_minus_one,solver_flag);
S=AdjustValueAtInfinityFromRattApp(S);
fprintf('Time to compute rat. app. %f\n',toc(t_rat_app));
if (plot_flag)
    PlotS(S);
    RaiseAndDraw;
end
 

