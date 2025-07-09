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
%line 391 <arl2.nw>
function rarl2Report = RARL2(data, N, maxIP, options)
t0 = cputime;
data = arl2Preprocess(data);
% Clear
arl2GetLMList;  
%line 399 <arl2.nw>
% Zero degree initialization
deg = 0;
[CSP0, sys0, J0] = arl2ZeroApprox(data);
arl2StoreLM(CSP0, sys0, J0, 0);
LM0 = arl2GetLMList; 
F0 = data.F0; 

nJ0 = sqrt(J0 / F0);
fprintf('-- Degree %d --  Local minimum: 1  J=%.7f\n', deg, nJ0);
%line 411 <arl2.nw>
% Increase degree upto N
lmlist = LM0; J = J0;
while deg<N
 LL = length(lmlist);
 deg = deg+1;
 if LL*maxIP>1, suffix='s'; else, suffix = ''; end
 fprintf('\n-- Degree %d --  Start with %d initial point%s\n', ...
 deg, LL*maxIP, suffix);
 printFlag = 0; 
 bestValue = Inf;
 for lmindex = 1:LL
  lm = lmlist{lmindex};
  ipList = ...
    arl2CreateIPList(lm, maxIP, data);    % Create initial points list
  for ipIndex = 1:length(ipList)
    ip = ipList{ipIndex};
    [sys, report, CSP] = ...
               arl2(deg, data, ip, options);  % Best approx starting from ip
    J = report.finalValue;
    printFlag = printFlag + 1; 
    if printFlag > 4, printFlag=1; fprintf('\n'); end
    nJ = sqrt(J / F0);
    fprintf('%14.7e ... ',nJ); 
    arl2StoreLM(CSP, sys, J, lmindex);        % Store if not yet found
    if nJ<bestValue, bestValue=nJ; end
  end
 end
 lmlist  = arl2GetLMList;                     % Prepare for next degree
 LM{deg} = lmlist;                            % Remember
 if length(lmlist)>1, suffix='a '; else, suffix = 'um'; end
 fprintf(...
 '\n-- Degree %d -- Found %d local minim%s -- Best relative error = %.7f\n', ...
    deg, length(lmlist), suffix, bestValue);
end
t1 = cputime;
rarl2Report.data     = data;
rarl2Report.n        = N;
rarl2Report.maxIP    = maxIP;
rarl2Report.data     = data;
rarl2Report.exectime = t1 - t0;
rarl2Report.LM0      = LM0;
rarl2Report.LM       = LM;
%line 456 <arl2.nw>
% ---------------- arl2StoreLM
function k = arl2StoreLM(CSP, sys, J, lmindex)
global arl2LMList
L = length(arl2LMList); 

for k = 1:L
 lm = arl2LMList{k};
 x = SSdistH2(sys, lm.sys)/sqrt(SSnH2(sys)*SSnH2(lm.sys));
 if arl2IsSmall(x, 1.e-3)
   lm.Father = [lm.Father k];
   arl2LMList{k} = lm;
   if J<lm.J, lm.J = J; end
   return
 end
end
k = L+1;
arl2LMList{k} = struct('CSP',CSP,'sys',sys,'J',J,'Father',lmindex);
%line 476 <arl2.nw>
% ---------------- arl2GetLMList
function lmlist  = arl2GetLMList
global arl2LMList
lmlist = arl2LMList; 
arl2LMList = {};
%line 484 <arl2.nw>
% ---------------- arl2ZeroApprox
function [CSP0, sys0, J0] = arl2ZeroApprox(data)
switch data.type
case 'sys'
 [a,b,c,d] = SSdata(data.value);
 [p,m] = size(d);
case 'sample'
 [p,m,LL]   = size(data.value);
case 'coef'
 [p,m,LL]   = size(data.value);
otherwise
 mesg = sprintf('arl2ZeroApprox: not implemented type %s', data.type);
 error(mesg)
end

p0 = min(p,m);
D0 = eye(p0,p0);
CSP0 = arl2CreateCSP(zeros(p0,0),[],zeros(p0,0),D0);
allpass0 = arl2CSP2Allpass(CSP0);
[J0, sys0] = arl2UserFunction(allpass0, data);
%line 507 <arl2.nw>
function CSP = arl2CreateCSP(u,w,y,D)
chart.u=u; chart.w=w; chart.D0 = eye(size(D));
CSP.chart=chart;
CSP.y=y;
CSP.D=D;
