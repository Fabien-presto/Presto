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
%line 309 <example.nw>
function [sys, report] = IRCOM(params)

close all
sys =[]; report = [];

% MacMillan degree
n = params.n;

% Data dimensions
m = params.m;
p = params.p;

% Data filenames (Fourier coefficients)
for lig = 1:p, for col = 1:m
  DataFiles{lig, col} = sprintf('%s%d%d%s',params.prefix, lig, col, params.postfix);
end, end

[coeff, values ] = ShowPrint(DataFiles);

m1 = sprintf('Fourier coefficients in %s', params.prefix );
m2 = 'Data: Nyquist and Bode locus (module)';
if ~ arl2Continue({m2, m1}),  return, end
%line 334 <example.nw>
%========================================================================
%                                ARL2 Part
%========================================================================

% Prepare ARL2 data struct 
data.type  = 'coef';
data.value = coeff;

% Random initial n-m-p system
%sys0 = arl2CreateStableSystem(n,m,p);
disp('Find initial system')
nbCoef = size(coeff,3);
nbCoef=min(nbCoef,50);
y = zeros(nbCoef,4);
for coefIndex = 1:nbCoef-1
  x = coeff(:,:,coefIndex);
  y(coefIndex+1,:) = x(:)';
end
[a,b,c,d,totbnd,hsv] = imp2ss(y, n, 2, 2);
sys0 = SScreate(a,b,c,d)


% Prepare figures
close all
% Data/approximation figure
global SSplotFigure
H = 100*p+50; W = 100*m+100;
pos = [ 10   10   W  H];
SSplotFigure = ...
 figure('Position',pos,'Name','Data & approximant','Visible','off');
SSplotFigure = [];
SSplotMat(values,'b');

SSplot(sys0,'g');
% m1 = 'Nyquist: data and initial system';
% if ~ arl2Continue({m1}),  return, end

% Minimizer options
options = { 			...
  'TolX',		eps,	...
  'TolFun',		eps,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'on'   ...
};

% Call minimizer
allpass0 = arl2AllpassFactor(sys0);
data.valinf = 0;

[sys, report] = arl2(n, data, allpass0, options);

if isempty(report), return, end

%if ~ arl2Continue('Refresh figure'), return, end
clf
%SSplotMat(values,'b');
%SSplot(sys,'r');
ShowPrint(DataFiles, sys);
if ~ arl2Continue('OK ?'), return, end

clc;
disp(report);

fprintf('\n');
Fplus = 1;
relError = report.finalValue  / Fplus;
fprintf('Initial value  : %.3e\n',Fplus);
fprintf('Final value    : %.3e\n',report.finalValue);
fprintf('Relative error : %.3e\n',relError);
fprintf('\n');
fprintf('Hessian (cond) : %.3e\n',cond(report.HESSIAN));


[A,B,C,D] = SSdata(sys);
fprintf('Poles approximant\n');
matPrint(sort(eig(A)));

save '../Data/cnes_may95r_sys8' sys

%line 416 <example.nw>
function [coeff, values ] = ShowPrint(DataFiles, sys)

close all
global SSplotFigure

% Data dimensions
[p, m] =size(DataFiles);

if nargin>1
  [a,b,c,d] = SSdata(sys);
  deg = size(a, 2);
  [P, M] = size(c*b);
  
  if any([P, M] ~= [p, m]), warning('Dimension mismatch'),
    M, P, m, p
sys
sys.d
    return 
  end
end

% Read data
for lig = 1:p, for col = 1:m
 file = DataFiles{lig,col};
 [RE , IM] = textread(file,'%f %f','commentstyle','shell');
 Data = RE + i * IM;
 % Split analytic / anti-analytic parts
 CoeffPlus (lig, col, :)  = Data;
 valInf(lig,col) = GetValInf(file);
end,end


% Plot data
N = length(Data);
coeff(:,:,1) = zeros(p,m);
coeff(:,:,2:N+1) = CoeffPlus;
values = fft(coeff,2*N+1,3);

% Prepare  figures

% ARL2 Data/approximation figure
if 0
H = 200*p+50; W = 200*m+100;
pos = [ 10   600   W  H];
figure('Position',pos,'Name','ARL2 data & approximant');
SSplotFigure = gcf;
SSplotMat(values,'b');
SSplot(sys,'r');
FigureTitle('ARL2 data and approximant - MacMillan degree: 8');
print -depsc2 ARL2.eps
end

% CNES figure (TO BE MODIFIED)
H = 200*p+50; W = 200*m+100;
pos = [ 600   600   W  H];
figure('Position',pos,'Name','CNES data & approximant');
SSplotFigure = gcf;

%% Values and sys should be modified

z = exp(2*pi*i*(0:2*N)/(2*N+1));
for lig = 1:p, for col = 1:m
cnesValue(lig,col,:) = ...
 valInf(lig,col) + fftshift((transpose(squeeze(values(lig,col,:))).*(z - 1)));
end, end

SSplotMat(cnesValue,'b');



if nargin > 1
  
  resp = SSresp(sys,2*N+1);
  z = exp(i*resp.theta);
  for lig = 1:p, for col = 1:m
  sysValue(lig,col,:) = ...
   valInf(lig,col) + (transpose(squeeze(resp.value(lig,col,:))).*(z - 1));
  end, end
  SSplotMat(sysValue,'r');
  
  dT = 2*pi/(2*N+1);
  errValue = cnesValue - sysValue;
  errNorm = sum(sum(sum(errValue.*conj(errValue)))) * dT;
  cnesNorm = sum(sum(sum(cnesValue.*conj(cnesValue)))) * dT;
  sysNorm = sum(sum(sum(sysValue.*conj(sysValue)))) * dT;
  
  relErr = errNorm/cnesNorm;
  
  titleFormat = '%s  -  CNES data and approximant  -  Degree: %d  -  Error: %.1e  (%.1f%%)';
  
  FigureTitle(sprintf(titleFormat, 'Nyquist', deg, relErr, 100*sqrt(relErr) ));
  print -depsc2 nyquistCNES.eps
end

H = 200*p+50; W = 200*m+100;
pos = [ 0   600   W  H];
figure('Position',pos,'Name','CNES data & approximant');
SSplotFigure = gcf;

SSbodeMat(cnesValue,'b');
if nargin > 1
  SSbodeMat(sysValue,'r');

FigureTitle(sprintf(titleFormat, 'Bode', deg, relErr, 100*sqrt(relErr) ));
print -depsc2 bodeCNES.eps
end



function FigureTitle(text)
fig = gcf;
pos = get(fig,'Position');
uicontrol(fig,...
	  'Style','Text',...
	  'String',text,...
	  'BackgroundColor','white',...
	  'Position',[0 0 pos(3) 20]...
	  );


function val = GetValInf(file)

fid = fopen(file);
if fid < 0, error([file , ' : cannot open']), end
while 1
 s = fgets(fid);
 pattern = '##/ re_inf =';
 if strncmp(s, pattern, length(pattern))
  re = sscanf(s, [pattern, '%f']);
 end
 pattern = '##/ im_inf =';
 if strncmp(s, pattern, length(pattern))
  im = sscanf(s, [pattern, '%f']);
  break
 end
end
 
val = re + i * im;

