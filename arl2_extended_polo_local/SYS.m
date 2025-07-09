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
%line 149 <example.nw>
function [sys, report, rp] = SYS(params,varargin)

clc; fprintf('\nCheck ARL2 minimization\n');
close all
%line 155 <example.nw>
global runParams
% Defaults
if isempty(runParams)
 runParams.n  = 4;
 runParams.m  = 2;
 runParams.p  = 3;
 runParams.n0 = 4;
end
% Overload ?
if nargin > 0
  if isfield(params,'n'),       runParams.n       = params.n  ; 	end  
  if isfield(params,'m'),       runParams.m       = params.m  ;		end  
  if isfield(params,'p'),       runParams.p       = params.p  ; 	end  
  if isfield(params,'n0'),      runParams.n0      = params.n0 ; 	end 
end
% Extract
n =   runParams.n  ;
m =   runParams.m  ;
p =   runParams.p  ;
n0 =  runParams.n0 ;
rp = runParams;

global SSplotFigure
H = 100*p+50; W = 100*m+100;
pos = [ 10   10   W  H];
SSplotFigure = figure('Position',pos,'Name','Data and approximation');
%line 183 <example.nw>
% Create unknown random system
TheSys = arl2CreateStableSystem(n,m,p);
SSplot(TheSys);
clc; SSdisplay(TheSys);
mesg = sprintf('"Unknown" system TheSys: n=%d,m=%d,p=%d\n', n, m, p);
if ~arl2Continue(mesg) , return, end
%line 191 <example.nw>
% Create initial random system
allpass0 = arl2CreateAllpass(n0,min(m,p));
data.type  = 'sys';
data.value = TheSys;
[err0, sys0]  = arl2UserFunction(allpass0, data);
SSplot(sys0,'g');
clc; SSdisplay(sys0);
mesg = sprintf('Initial random system Sys0: n=%d,m=%d,p=%d\n', n0, m, p);
if ~arl2Continue(mesg) , return, end
%line 202 <example.nw>
% Norms and distance
F0 = SSnH2(TheSys);        
m1 = sprintf('Unknown system TheSys H2 norm: %f\n',sqrt(F0)); 

F1 = SSnH2(sys0)  ;        
m2 = sprintf('Initial system Sys0   H2 norm: %f\n',sqrt(F1));

d0=SSdistH2(TheSys, sys0); 
m3 = sprintf('Systems H2 distance: %f (FVAL=%f)\n', sqrt(d0), d0);

arl2Assert(arl2IsSmall(abs(d0-err0)));
%line 215 <example.nw>
mesg = sprintf('\nFirst test: system exact matching\n');
if arl2Continue({m1,m2,m3,mesg})
%if 0
 % Try to approximate

 relTol = 1.e-12;
 absTol = relTol * F0;

 % Minimizer options
 options = { 			...
   'TolFun',		absTol,	...
   'GradObj',  		'on' , 	...
   'DerivativeCheck',	'on'   ...
 };
 [sys, report] = arl2(n, data, allpass0, varargin);
 
 SSplot(sys,'r');
 if ~isempty(sys) & ~isempty(report) 
  fprintf('\nFinal distance: %e\n\n', SSdistH2(TheSys, sys)) ;
  report.relTol       = relTol;
  report.absTol       = absTol;
  report.initialValue = F0;
  report.n = n;
  report.m = m;
  report.p = p;
  registerReport(report);
 end
end
%line 245 <example.nw>
mesg = sprintf('\nSecond test: sampled data matching\n');

if arl2Continue(mesg)
 data = SSresp(TheSys,50);
 data.type = 'sample';
 figure(SSplotFigure); clf
 SSplot(TheSys,'b');
 SSplotMat(data.value, 'k.');
 SSplot(sys0,'g');

 % Try to approximate
 
 newData = arl2Preprocess(data);
 fprintf('\n TheSys and noisy data matching ? %e %e\n', ...
 abs(newData.F0 - F0), norm(newData.D - TheSys.d));

 
 
 relTol = 1.e-12;
 absTol = relTol * F0;

 % Minimizer options
 options = { 			...
   'TolFun',		absTol,	...
   'GradObj',  		'on' , 	...
   'DerivativeCheck',	'on'   ...
 };
 [sys, report] = arl2(n, data, allpass0, options);
 SSplot(sys,'r');

 if ~isempty(sys) & ~isempty(report) 
  fprintf('\nFinal distance: %e\n', SSdistH2(TheSys, sys)) ;
  report.relTol       = relTol;
  report.absTol       = absTol;
  report.initialValue = F0;
  report.n = n;
  report.m = m;
  report.p = p;
  registerReport(report);
 end
end
%line 288 <example.nw>
function registerReport(r)

if isempty(r), return, end

month = str2num(datestr(now,'mm'));
day   = str2num(datestr(now,'dd'));
fid = fopen('arl2LOGFILE','a');
p = min(r.p, r.m);
problemSize = (2*r.n + p)*p;
formatstr = ...
 '%02d/%02d %2d %2d %2d %3d %8.1e %8.1e %8.1e %8.1e %5d %5d %7.1f\n';
fprintf(fid, formatstr, ...
 day, month, ...
 r.n, r.m, r.p, problemSize, ...
 r.initialValue,  r.finalValue, ...
 r.relTol, r.absTol, ...
 r.iterations, r.funcCount, r.execTime);
fclose(fid);
