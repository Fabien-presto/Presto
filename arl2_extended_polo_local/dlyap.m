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
function x = dlyap(a,b,c)
%DLYAP	Discrete Lyapunov equation solver.
%
%    X = DLYAP(A,Q) solves the discrete Lyapunov equation:
%
%		A*X*A' - X + Q = 0
%
%    See also  LYAP.

%	J.N. Little 2-1-86, AFP 7-28-94
%	Copyright (c) 1986-1999 The Mathworks, Inc. All Rights Reserved.
%	$Revision: 1.2 $  $Date: 2001/01/17 16:47:31 $

% How to prove the following conversion is true.  Re: show that if
%         (1) Ad X Ad' + Cd = X             Discrete lyaponuv eqn
%         (2) Ac = inv(Ad + I) (Ad - I)     From dlyap
%         (3) Cc = (I - Ac) Cd (I - Ac')/2  From dlyap
% Then
%         (4) Ac X + X Ac' + Cc = 0         Continuous lyapunov
% 
% Step 1) Substitute (2) into (3)
%         Use identity 2*inv(M+I) = I - inv(M+I)*(M-I) 
%                                 = I - (M-I)*inv(M-I) to show
%         (5) Cc = 4*inv(Ad + I)*Cd*inv(Ad' + I)
% Step 2) Substitute (2) and (5) into (4)
% Step 3) Replace (Ad - I) with (Ad + I -2I)
%         Replace (Ad' - I) with (Ad' + I -2I)
% Step 4) Multiply through and simplify to get
%         X -inv(Ad+I)*X -X*inv(Ad'+I) +inv(Ad+I)*Cd*inv(Ad'+I) = 0
% Step 5) Left multiply by (Ad + I) and right multiply by (Ad' + I)
% Step 6) Simplify to (1)

ni = nargin;

if ni==2,
   c = b;
   b = a';
end

[m,n] = size(b);

a = (a+eye(size(a)))\(a-eye(size(a)));
b = (b+eye(size(b)))\(b-eye(size(b)));
c = (eye(size(a))-a)*c*(eye(size(b))-b)/2;
x = lyap(a,b,c);

% end dlyap
