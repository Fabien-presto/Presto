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
%line 558 <example.nw>
function report =  Bench(maxIP)
%line 562 <example.nw>
% RARL2 directory
ARL2DIR = '/0/home/jpm/Work/Arl2/Noweb';
addpath(ARL2DIR);

% Initial points per local minimum: default 4
if nargin<1, maxIP = 4; end

mesg = 'RARL2 -- Choose reference system';
System = questdlg(...
     mesg,'RARL2 system approximation','HB61', 'GloverExample', 'ex2to1', 'HB61' );

% Load system (discrete-time) and dimensions
TheSys  = feval(System);
[n,m,p] = SSdim(TheSys);
nS = sqrt(SSnH2(TheSys));
TheSys.c = TheSys.c / nS;
TheSys.d = TheSys.d / nS;
%line 582 <example.nw>
clc; close all; 
fprintf('\nCheck ARL2 minimization - System: %s - H2-norm = %e\n\n', ...
 System, sqrt(SSnH2(TheSys)));

% Plot and display discrete-time system, ask to continue ?
ShowSys(TheSys);
mesg = sprintf('Reference system (%s): n=%d  m=%d  p=%d\n', System, n, m, p);
if ~arl2Continue(mesg) , return, end
%line 593 <example.nw>
% Prepare data
data.type  = 'sys';
data.value = TheSys;

% Tolerances and minimizer options
options = { 			...
  'Display',            'off',  ...
  'TolFun',		eps,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'off'    ...
};
%line 607 <example.nw>
% Here we go
report = RARL2(data, n, maxIP, options )
%line 612 <example.nw>
% Reporting
report.message  = ['RARL2 - System:  ' System];
Title = 'Save report ?';
FileSpecs = sprintf('REPORTS/%s.mat',System);
[FileName, PathName] = uiputfile(FileSpecs, Title);
if FileName ~= 0,  save( [PathName FileName],'report' ); end

ShowReport(report)

%line 624 <example.nw>
% ===================================================================== 
% ======================= Utilities and systems =======================
% ===================================================================== 
function ShowSys(sys);
global SSplotFigure
[n,m,p] = SSdim(sys);
H = 100*p+50; W = 100*m+100;
pos = [ 10   10   W  H];
SSplotFigure = figure('Position',pos,'Name','Data and approximation');
SSplot(sys);
SSdisplay(sys);
%line 638 <example.nw>
% ===================================================================== GloverSystem
function [ DiscreteSystem, ContinuousSystem ] = GloverExample

% Continuous system
a = [
 0.00  1.00  0.00   0.00   0.00  0.00  0.00  0.00    0.00    0.00    0.00   0.00
-0.20 -1.15  0.00   0.00   0.00  0.00  0.00  0.00    0.00    0.00    0.00   0.00
 0.00  0.00  0.00   1.00   0.00  0.00  0.00  0.00    0.00    0.00    0.00   0.00
 0.00  0.00  0.00   0.00   1.00  0.00  0.00  0.00    0.00    0.00    0.00   0.00
 0.00  0.00 -2.36 -13.60 -12.80  0.00  0.00  0.00    0.00    0.00    0.00   0.00
 0.00  0.00  0.00   0.00   0.00  0.00  1.00  0.00    0.00    0.00    0.00   0.00
 0.00  0.00  0.00   0.00   0.00  0.00  0.00  1.00    0.00    0.00    0.00   0.00
 0.00  0.00  0.00   0.00   0.00 -1.62 -9.40 -9.15    0.00    0.00    0.00   0.00
 0.00  0.00  0.00   0.00   0.00  0.00  0.00  0.00    0.00    1.00    0.00   0.00
 0.00  0.00  0.00   0.00   0.00  0.00  0.00  0.00    0.00    0.00    1.00   0.00
 0.00  0.00  0.00   0.00   0.00  0.00  0.00  0.00    0.00    0.00    0.00   1.00
 0.00  0.00  0.00   0.00   0.00  0.00  0.00  0.00 -188.00 -111.60 -116.40 -20.80
];

b = [
   0.00    0.00
   1.04    4.15
   0.00    0.00
   0.00    0.00
  -1.79    2.68
   0.00    0.00
   0.00    0.00
   1.04    4.15
   0.00    0.00
   0.00    0.00
   0.00    0.00
  -1.79    2.68
];

c = [
  0.26  0.81  -1.42  -15.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00
  0.00  0.00   0.00    0.00  0.00  4.90  2.12  1.95  9.35 25.80  7.14  0.00
];

d = [
 0.00    0.00
 0.00    0.00
];

% Discrete system
I  = eye(size(a));
AI = inv(I-a);
A  = -(I+a)*AI;
B  = b;
C  = c*AI;
D  = d;

ContinuousSystem = SScreate(a,b,c,d);
DiscreteSystem   = SScreate(A,B,C,D);
%line 695 <example.nw>
% ===================================================================== HB61System
function [ DiscreteSystem, ContinuousSystem ] = HB61

% Continuous system
a = [
	0	0	0	-150
	1	0	0	-245
	0	1	0	-113
	0	0	1	-19
];

b = [
	4
	1
	0
	0
];

c = [	0	0	0	1];

d = 0;

% Discrete system
I=eye(size(a));
AI=inv(I-a);
A=-(I+a)*AI;
B=b;
C=c*AI;
D=d;

ContinuousSystem = SScreate(a,b,c,d);
DiscreteSystem   = SScreate(A,B,C,D);
%line 730 <example.nw>
% ===================================================================== ex2to1System
function DiscreteSystem  =  ex2to1

% Discrete system
m=1;
p=1;

% Systeme reel qui admet un approx complexe -2<b1<0 
a = [
 0	1	0
 0	0	1
 0	0	0
];

b = [
   -1
    0
    1
];

c = [1  0  0];

d = 0;


DiscreteSystem   = SScreate(a,b,c,d);


