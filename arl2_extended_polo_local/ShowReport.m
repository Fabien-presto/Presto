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
%line 761 <example.nw>
function varargout = ShowReport(varargin)

if nargin == 0 | ~ischar(varargin{1})                   % LAUNCH GUI
%line 767 <example.nw>
        if nargin == 0  , 
	    FileSpecs = 'REPORTS/*.mat'; Title = 'Choose report';
            [FileName, PathName] = uigetfile(FileSpecs, Title);
	    if FileName == 0, return, end;
	    load([PathName, FileName]);
        else,             report = varargin{1};
        end

	disp(report);

%line 780 <example.nw>
	close all
	warning off; 
	fig = openfig(mfilename);

	set(fig, 'Name', report.message);
	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
        handles.fig = fig;
        handles.report = report;
        guidata(fig, handles);

        DrawRefSys(handles, report);

	warning on;

%line 801 <example.nw>
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end
%line 817 <example.nw>
function DrawRefSys(handles, report)
htitle = handles.title;
hpopup = handles.popup;

maxDeg = length(report.LM);

degList = {'None'};
for deg = 1:maxDeg
degList = {degList{:}, sprintf('%d', deg)};
end
%degList = {degList{:}, 'All'};

set(hpopup,'String',degList);
ShowData(handles.fig, report.data);
set(htitle, 'String', sprintf('Reference system\n'));
%line 835 <example.nw>
function ShowReportCallback
handles = guidata(gcf);
num=get(handles.popup,'Value');
str = get(handles.popup,'String');
sel = str{num};
report = handles.report;

switch sel
case 'None'
 cla
 set(handles.title, 'String', sprintf('Reference system\n'));
 ShowData(handles.fig, report.data);
 return
otherwise
 k = str2num(sel);
      LL = length(report.LM{k});
      if LL>1, suffix='a '; else, suffix = 'um'; end
      bestJ=approxAndSys(handles, report, k);
      mesg = sprintf('Degree %2d - %2d minim%s - Best J=%.6f\n',k,LL,suffix,bestJ);
      set(handles.title, 'String', mesg)
end
%line 859 <example.nw>
function   bestJ = approxAndSys(handles, report, k)
 colors = {'r', 'g', 'k', 'y', 'm'};
  ShowData(handles.fig, report.data); 
  LL = length(report.LM{k});
  bestJ=Inf;
  for l = 1:length(report.LM{k})
  lm = report.LM{k}{l};
  if lm.J<bestJ, bestJ=lm.J; end
  col = colors{mod(l-1,length(colors))+1};
  ShowData(handles.fig, struct('type','sys','value',lm.sys), col, 'hold')
  end
  bestJ = sqrt(bestJ / report.data.F0 );
%line 874 <example.nw>
function ShowData(fig, data, varargin)
figure(fig)

if length(varargin)>=1, color=varargin{1}; else, color = 'b'; end
if length(varargin)>=2, holdFlag=varargin{2}; else, holdFlag='no'; end
%line 882 <example.nw>
switch data.type
case 'sys'
 sys = data.value;
 [n,m,p] = SSdim(sys);
 resps = SSresp(sys, 1000);
case 'coef'
 [p,m] = size(data.D);
 N = length(data.value); [p,m] = size(data.D);
 coeff(:,:,1) = data.D;  coeff(:,:,2:N+1) = data.value;
 resps.value = fft(coeff,2*N+2,3);
 resps.theta = 2*pi*(0:2*N+1)/(2*N+2);
case 'sample'
 [p,m] = size(data.D);
 resps = data;
end
%line 900 <example.nw>
 index = 0;
 for lig=1:p
 for col=1:m
  index = index + 1;
  subplot(p,m,index);
  if strcmp(holdFlag,'no'), cla, end
  resp=squeeze(resps.value(lig,col,:));
  plot(real(resp), imag(resp),color);
  axis equal
  hold on
 end
 end
 drawnow
