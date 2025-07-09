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
%line 918 <example.nw>
function RARL2Launcher

FileSpecs = 'DATA/*.mat'; Title = 'Select data';
[FileName, PathName] = uigetfile(FileSpecs, Title);
if FileName == 0, return, end
load([PathName, FileName]);
if isfield(data, 'name'), System = data.name; 
else, System = 'unknown'; end
%line 929 <example.nw>
fprintf('Processing %s%s\nName: %s\nData type: %s\n\n',...
	PathName, FileName, System, data.type);
Prompt = {'Approximant degree ?', 'Initial points number'};
Title  = ['Data: ' data.name];
LineNo = [1 40;1 40];
DefAns = {'8','4'};
Answer = inputdlg(Prompt,Title,LineNo,DefAns);
if isempty(Answer), return, end

n     = str2num(Answer{1});
maxIP = str2num(Answer{2});

% Tolerances and minimizer options
options = { 			...
  'Display',            'off',  ...
  'TolFun',		eps,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'off'    ...
};
%line 951 <example.nw>
% Here we go
fprintf('Mac-Millan degree: %d\n',  n);
fprintf('Starting points per local minimum at each degree: %d\n\n', maxIP);

report = RARL2( data, n, maxIP, options );

% Add message
report.message  = ['RARL2 approximation - Data: ' PathName FileName];

Title = 'Save report ?';
FileSpecs = sprintf('REPORTS/%s.mat',System);
[FileName, PathName] = uiputfile(FileSpecs, Title);
if FileName ~= 0,  save( [PathName FileName], 'report' ); end

ShowReport(report);

