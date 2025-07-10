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
% Get current directory
DemoDir=GetDir;

% Set presto in path
addpath([DemoDir,'presto']);

% Init presto (sets a bunch of paths)
PrestoInit;

% Loads the data (the data are after that in a structure called eight_cavity_dual_mode_filter)
eight_cavity_dual_mode_filter=ParseS('eight_cavity_dual_mode_filter.txt');

% Plot the nyquist plot of data
PlotS(eight_cavity_dual_mode_filter);
RaiseAndDraw;

% Plot the bode plot of data
PlotS(eight_cavity_dual_mode_filter,'b');
RaiseAndDraw;

% Indicates in the eight_cavity_dual_mode_filter structure where to take special parameter
% value. This is only necessary if one wants to use different parameters
% as the default ones.
eight_cavity_dual_mode_filter.params=@Params_for_eight_cavity_dual_mode_filter

% Compensate delays and freq shift - compute completion at infinity
Sr=CompensateDelayAndFreqShift(eight_cavity_dual_mode_filter,8);

% Plot the result
PlotS(Sr);
RaiseAndDraw;

% Compute Fourier coeffs and evaluate importance of non-causal part
Sr=CompFourierCoeffs(Sr);
PlotS(Sr);
RaiseAndDraw;

% Compute a rational approximant of the compensated and completed matrix 
% of MacMillan degree 8 - the last arguments indicates tha twe pass to
% the circle with an H^2 isometry
Sr=RatApp(Sr,8,1)
PlotS(Sr);
RaiseAndDraw;

% All those steps can be called by one command by: Sr=AllSteps(eight_cavity_dual_mode_filter,8)


% Compute the extended coupling matrix, the classical coupiling matrix, and the input/output loads of Sr
[M,Mred,R1,R2]=S2M(Sr);

% Compute the arrow form

disp('Forme en fleche: partie reelle');
real(M)
disp('Forme en fleche: partie imaginaire');
imag(M)

 
 

