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
function X = lyap(A, B, C)
%LYAP  Solve continuous-time Lyapunov equations.
%
%   X = LYAP(A,C) solves the special form of the Lyapunov matrix 
%   equation:
%
%           A*X + X*A' = -C
%
%   X = LYAP(A,B,C) solves the general form of the Lyapunov matrix
%   equation (also called Sylvester equation):
%
%           A*X + X*B = -C
%
%   See also  DLYAP.

%	S.N. Bangert 1-10-86
%	Copyright (c) 1986-1999 The Mathworks, Inc. All Rights Reserved.
%	$Revision: 1.1 $  $Date: 2000/09/07 06:56:01 $
%	Last revised JNL 3-24-88, AFP 9-3-95

ni = nargin;

if ni==2,
   C = B;
   B = A';
end

[ma,na] = size(A);
[mb,nb] = size(B);
[mc,nc] = size(C);

% A and B must be square and C must have the rows of A and columns of B
if (ma ~= na) | (mb ~= nb) | (mc ~= ma) | (nc ~= mb)
   error('Dimensions do not agree.')
elseif ma==0,
   X = [];
   return
end

% Perform schur decomposition on A (and convert to complex form)
[ua,ta] = schur(A);
[ua,ta] = rsf2csf(ua,ta);
if ni==2,
   % Schur decomposition of A' can be calculated from that of A.
   j = ma:-1:1;
   ub = ua(:,j);
   tb = ta(j,j)';
else
   % Perform schur decomposition on B (and convert to complex form)
   [ub,tb] = schur(B);
   [ub,tb] = rsf2csf(ub,tb);
end

% Check all combinations of ta(i,i)+tb(j,j) for zero
p1 = diag(ta).'; % Use .' instead of ' in case A and B are not real
p1 = p1(ones(mb,1),:);
p2 = diag(tb);
p2 = p2(:,ones(ma,1));
sum = abs(p1) + abs(p2);
if any(any(sum == 0)) | any(any(abs(p1 + p2) < 1000*eps*sum))
   error('Solution does not exist or is not unique.');
end

% Transform C
ucu = -ua'*C*ub;

% Solve for first column of transformed solution
y = zeros(ma,mb);
ema = eye(ma);
y(:,1) = (ta+ema*tb(1,1))\ucu(:,1);

% Solve for remaining columns of transformed solution
for k=2:mb,
   km1 = 1:(k-1);
   y(:,k) = (ta+ema*tb(k,k))\(ucu(:,k)-y(:,km1)*tb(km1,k));
end

% Find untransformed solution 
X = ua*y*ub';

% Ignore complex part if real inputs (better be small)
if isreal(A) & isreal(B) & isreal(C),
   X = real(X);
end

% Force X to be symmetric if ni==2 and C is symmetric
if (ni==2) & isequal(C,C'),
   X = (X+X')/2;
end

% end lyap
