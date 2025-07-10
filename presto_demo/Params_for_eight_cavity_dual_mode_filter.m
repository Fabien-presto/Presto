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
function DVC=Params_for_eight_cavity_dual_mode_filter()
% Returns a parameter structure containing Default values and Constants values for all Presto methods 
% For a global change of those values you can edit the file PrestoDefaultValuesAndConstants.m
% For a local overloading of those values define a 'params' field for your S structure
% pointing to a function handle. This handle should return (after evaluation) a parameter structure 
% with the parameters you want to overload.
% ex: S.params=@GetMyDefaultValuesAndConstants
%     AllSteps(S)
% In this example AllSteps will take into acount the overloaded values defined in GetMyDefaultValuesAndConstants.m



%******************************************************
% DVC for CompensateDelayAndFreqShift
%******************************************************

%     ************  Completion and delay determination parameters *******

DVC.CDAFS.delay_range=[-0.1:.001:0.1];
% Fixes the range of different delays used in the delay determination procedure
% Recall that delay determination is done by exhaustive trial, i.e each point in 
% delay_range leeds to an evaluation

DVC.CDAFS.poly_order=4;
% Fixes the degree of the serie (in z) used to compute completion at infinity
% The assumption is that this number can be taken small, as the filter behaviour is smooth 
% at infinity

DVC.CDAFS.omega_lim=2.5;
% For frequency greater than omega_lim and smaller thant -omega.lim the algorithm tries to fit a serie 
% expansion in z to the data for the diagonal terms (S11, S22). This is how the delay determination 
% works. So for w>omega.lim the data should behave quite smoothly, allowing delay determination using 
% measurements in this frequency range.
% Important: note that omega.lim should be choosen so that there are still enough measuremnts points 
% for |w|>omega.lim to determine a meaningfull serie expansion in 1/s.

DVC.CDAFS.causal_bound(1,1)=0.5;
% For the (1,1) and (2,2) element a constrained probleme is solved to determine best completion. In this 
% problem the distance to causal functions is minimized so as to be less than the causal_bound if possible.
% This bound is expressed in %

DVC.CDAFS.causal_bound(2,2)=0.5;
% Same as prec.

DVC.CDAFS.modulus_factor(1,1)=1.01;
% The modulus of the completion is controled so as to be less than: modulus_factor * max(abs(S(1,1))  

DVC.CDAFS.modulus_factor(2,2)=1.01;
% Same as prec.







 

