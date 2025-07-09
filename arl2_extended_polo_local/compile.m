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
%line 407 <utilities.nw>
function compile(program, args)

if isunix
 [status, result] = unix('uname'); 
 System=result(1:5); 
else
 System = 'DOS'; 
end

switch System

 case 'SunOS'
   Mfile = sprintf('%s.m',program);   
   Mexfile = sprintf('%s.exe',program);   
   fprintf('Compiling %s to %s ... ', Mfile, Mexfile);
   command = sprintf('mcc -m %s -o %s',program,Mexfile);
   eval(command);
   fprintf('Done\n');
   command = sprintf('%s %s',Mexfile,args);
   fprintf('\nLaunch executable %s\n',command);
   unix(command);
 
 case 'Linux'
   disp('Sorry, problems with libc6 under Linux');
 
 otherwise
   disp(['Unknown system : ' , System ]);

end
