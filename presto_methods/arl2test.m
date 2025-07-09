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
function [S]=RatApp(S,n); 
addpath /net/miaou/lib/arl2/Matlab ;
addpath ../arl2polo;
assert(isfield(S,'fourier'),'Compute fourier coeffs first');
for k=1:2
    for j=1:2
        CoeffPlus(k,j,:)=conj(S.fourier{k,j}(2:end));        
    end
end
data.value = CoeffPlus;
data.type  = 'coef';
sys0 = arl2CreateStableSystem(n,2,2);
options = { 			...
  'TolFun',		eps,	...
  'GradObj',  		'on' , 	...
  'DerivativeCheck',	'on'   ...
};
% Schur parameters monitoring 
arl2InitMonitor([], 10, 40, 4); 
[sys, report] = arl2(10, data, sys0, options);
close all;
S.rat_sys=sys;
%resp = SSresp(sys,1000);
%theta=linspace(0,2*pi,100);
%v=conj((freqresp(ss(sys.a,sys.b,sys.c,sys.d,1),theta)))./exp(i*transpose(theta));