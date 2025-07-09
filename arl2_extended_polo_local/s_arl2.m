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
%line 8 <arl2.nw>
 
%line 12 <arl2.nw>
function [sys, report, newCSP] = arl2(finalDegree, userData, sys0 , varargin)

sys    = [];
report = [];

data    = arl2Preprocess(userData);
options = arl2CustomizeOptions(varargin{:});
if isfield(userData , 'F0')
 options.TolFun = options.TolFun * userData.F0;
end

%mesg = sprintf('\nCurrent options:\n');
%clc; disp(options);
%if ~arl2Continue(mesg) , return, end
%clc

%line 56 <arl2.nw>
% Always compute the allpas factor - Bug corrected with respect to previous
% code
% Encode initial system
%if arl2IsAllpass(sys0)
%    allpass = sys0;
% else
allpass = arl2AllpassFactor(sys0);
%end
CSP0    = arl2Allpass2CSP(allpass);

%SSdisplay(allpass);
%arl2CSPdisplay(CSP0);
%line 70 <arl2.nw>
% Minimize
T = cputime;
[newCSP, OUTPUT, HESSIAN, Xad, sys] = arl2MinimizePsi(finalDegree, options, CSP0, data);
OUTPUT.execTime = cputime-T;
OUTPUT.HESSIAN = HESSIAN;  
%line 378 <arl2.nw>
% Last build output arguments: system and report
allpass = arl2CSP2Allpass(newCSP);
report = OUTPUT;
report.Xad = Xad;
%line 9 <arl2.nw>
 
%line 30 <arl2.nw>
function options = arl2CustomizeOptions(optionList)

nbOptions = length(optionList);
if mod(nbOptions,2) ~= 0, nbOptions, optionList, error('Usage: option, value'); end

% Default options
options = optimset(...
		'MaxIter'	 	, 100		, ...
	        'Display'	 	, 'off'	        , ...
		'TolFun'	 	, eps		, ...
		'GradObj'	 	, 'on'		, ...
		'GradConstr'	 	, 'on'		, ...
		'DerivativeCheck'	, 'off'     ...     
	);
for index = 1:2:nbOptions
 currOption = optionList{index};
 currValue  = optionList{index+1};
 if isfield(options, currOption)
  options = optimset(options, currOption, currValue);
 else
  mesg = sprintf('%s: bad option', currOption );  error(mesg);
 end
end
