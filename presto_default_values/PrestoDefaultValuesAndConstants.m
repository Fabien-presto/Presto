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
function DVC=PrestoDefaultValuesAndConstants()
% DVC=PrestoDefaultValuesAndConstants()
% Returns a parameter structure containing Default values and Constants values for all Presto methods 
% For a global change of those values you can edit the file PrestoDefaultValuesAndConstants.m
% For a local overloading of those values define a 'params' field for your S structure
% pointing to a function handle. This handle should return (after evaluation) a parameter structure 
% with the parameters you want to overload.
% ex: S.params=@GetMyDefaultValuesAndConstants
%     AllSteps(S)
% In this example AllSteps will take into acount the overloaded values defined in GetMyDefaultValuesAndConstants.m



%***************************************************************
% DVC for AllSteps
%***************************************************************

%Constants

% Default values
DVC.AS.plot_flag=1; 
% inficates if results have to be ploted after each steps 
% poss. values: {0,1}

DVC.AS.solver_flag='RARL2'; 
% fixes which rational approximation engine is used by default
% poss. values: 'MRARL' 'RARL2' and 'HARL2'
% Important: HARL2 runs a binary file which has to be compatible with your system

DVC.AS.zeros_at_inf=2;  
% Specifies the annulation order at infinity for the diagonal terms S12 and S21
% Is used when computing completion of the data 
           
%******************************************************
% DVC for CompensateDelayAndFreqShift
%******************************************************

%Constants and magic numbers 

DVC.CDAFS.number_of_fourier_coeffs=250;
% The number of fourier coeffs used in contrained completion computations

DVC.CDAFS.number_of_control_points=20;
% The number of control points used to control modulus of completion

DVC.CDAFS.error_lim=0.05;
% Maximal error observed when trying to fit serie expansion to data. If the error is greater a warning 
% message is released.

DVC.CDAFS.sign_in_zero_trans_12=1; 
% Normalisation used to fix the sign of the imaginary part of S12 at zero frequency 
% This is used for normalisation reasons as the frequency shift is defined up to Pi
% If you think the Nyquist plot of S12 and S21 should be the other way arround change this
% value to its opposit.
% Poss. values: {-1,1}

DVC.CDAFS.iso_flag=0;
% This flag specifies the isometry used to estimate the non-causal part
% in the constrained problems solved for (1,1) and (2,2). 
% Poss. values: {0,1} 0: H^{\infty} isom, 1: H^2 isometry  


% ***   Default values ***
DVC.CDAFS.strategy='causal';
% This specifies which strategy is used to determine a delay compensation 
% Poss: {'causal','border_match'}

DVC.CDAFS.plot_flag=0;
% Poss. values: {0,1}

DVC.CDAFS.zeros_at_inf=2;
% See AllSteps def. values.




%     ************  Completion and delay determination other constants *******

DVC.CDAFS.delay_range=[-0.1:.001:0.1];

% Fixes the range of different delays used in the delay determination procedure
% Recall that delay determination is done by exhaustive trial, i.e each point in 
% delay_range leeds to an evaluation

DVC.CDAFS.poly_order=3;
% Fixes the degree of the serie (in z) used to compute completion at infinity
% The assumption is that this number can be taken small, as the filter behaviour is smooth 
% at infinity
% poly order can also be defined as an array, ex.
%DVC.CDAFS.poly_order(1,1)=10;
%DVC.CDAFS.poly_order(2,2)=10;
%DVC.CDAFS.poly_order(1,2)=3;
%DVC.CDAFS.poly_order(2,1)=3;
% In this case all four entries must be defined

DVC.CDAFS.omega_lim=2.5;
% For frequency greater than omega_lim and smaller thant -omega.lim the algorithm tries to fit a serie 
% expansion in z to the data for the diagonal terms (S11, S22). This is how the delay determination 
% works. So for w>omega.lim the data should behave quite smoothly, allowing delay determination using 
% measurements in this frequency range.
% Important: note that omega.lim should be choosen so that there are still enough measuremnts points 
% for |w|>omega.lim to determine a meaningfull serie expansion in z.

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



%******************************************************
% DVC for CompFourierCoeffs
%******************************************************

% Constants

DVC.CFC.n_fourier=250;
% Number of positive indexed fourier coeff. computed (i.e the algorithm computes indeed 2*n_fourier-1 coeffs)

DVC.CFC.n_comp=1000; 
% Number of discretisation points taken for evaluation of the completion (outside measurements range) to compute
% four. coeffs.

DVC.CFC.n_spline=2000;
% Number of discretisation points taken for evaluationof splined measurements to compute four. coeffs.

%******************************************************
% DVC for RatApp
%******************************************************

% Default values



DVC.RA.solver_flag='RARL2'; 
% fixes which rational approximation engine is used by default
% poss. values: 'MRARL2' 'RARL2' and 'HARL2' 
% Important: HARL2 runs a binary file which has to be compatible with your system

DVC.RA.divide_by_z_minus_one=0; 
% fixes the type of isometry used to pass to the circle - 
% poss. values: {0,1} - 0 is L{\infty} isom. 1 is is L^2 isom. 

DVC.RA.shift_adjustment=1;
% Additional correction of the fix shift for S12 and S21 during the rational approximation phase

DVC.RA.post_optimization=1;
% Start post optimisation with a l2 pointwise criterium at the end of classical rational approximation
% Can only be used with solver RARL2
% poss. values: {0,1}

DVC.RA.shift_adjustment_post=1;
% Additional correction of the fix shift for S12 and S21 during the post_rational approximation phase

% Stops ptimization if variation in the objevtive value is less than
% treshold - in the core optimization step
DVC.RA.func_stop_crit_optim=1e-7;

% Stops ptimization if variation in the objevtive value is less than
% treshold - in the post optimization step
DVC.RA.func_stop_crit_post_optim=1e-7;

%******************************************************
% DVC for PlotS
%******************************************************

% Constants

DVC.PS.max_w_for_nyq=500;
% Indicates for nyquist plot the maximal freq. to be ploted

DVC.PS.supp_ratio_for_bode=1.5;
% Indicates for bode plot the maximal ratio between ploted range and measurements ranges

% Default values

DVC.PS.legend_flag=1;
% Indicates if a legend is to be ploted for each curve (data, completion etc..)
% poss. values: {0,1}

DVC.PS.graphic_type_flag='n';
% Indicates if a Nyquist or a Bode plot is expected 
% poss. values: {'n','b'}



 

