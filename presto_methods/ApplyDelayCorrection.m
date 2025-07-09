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
function [Sout]=ApplyDelayCorrection(Sin,Sref)
% [Sout]=ApplyDelayCorrection(Sin,Sref)
% Applies the phase correction Sref.a and Sref.b to the data Sin

a=Sref.a;
b=Sref.b;
Sout=Sin;
for k=1:2
    for j=1:2
        Sout.value(k,j,:)=squeeze(Sout.value(k,j,:)).*(exp(a(k,j)*1i*Sout.freq));
        Sout.value(k,j,:)=Sout.value(k,j,:)*exp(1i*b(k,j)); 
    end
end
if (~isfield(Sin,'a'))
    Sout.a=zeros(2,2);
end
if (~isfield(Sin,'b'))
    Sout.b=zeros(2,2);
end
Sout.a=Sout.a+a;
Sout.b=Sout.b+b;

 

