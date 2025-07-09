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
function [sys,hsv]  = fc2ir(fc, nx)
% [sys, totbnd]  = fc2ir(fc, nx)
% 
% Convert Fourier coefficients to internal realization using AAK-theory
%

% Find out system dimensions
switch ndims(fc)
case 3
 [p, m, nc] = size(fc);
case 2
 [p , nc]  = size(fc);
 if min([p,nc])~=1
   mesg = sprintf('[%d %d]: ambiguous Fourier coefficients  dimensions', [p, nc]); 
   error(mesg);
 end
 fc = reshape(fc(:), [1,1, length(fc)]);
 [p, m, nc] = size(fc);
otherwise
 mesg = sprintf('Fourier coefficients: too many dimensions (%d)', ndims(fc)); 
 error(mesg); 
end

% Get D matrix D=H(0)
d = 0;
y=[];
for k = 1:nc
 y = [y fc(:,:,k)];
end

%disp('% Build Hankel matrix ')
if m==1 & p==1,
  n2 = ceil(length(y)/2);
  h = hankel(y(1:n2),y(n2:end));
else
   h = []; 
   n = floor((nc+1)/2);
   for k=1:n
      index = (k-1)*m+1 : (n-1+k)*m; 
      h = [h;y(:,index)];
   end
end
[rh,ch]=size(h);

%disp('% Compute the SVD of the Hankel')
[u,hsv,v] = svd(h);
hsv=diag(hsv);
ns=length(hsv);

tail=[conv(ones(ns,1)',hsv') 0];
tail=tail(ns:ns+ns);

if nargin < 2
%disp('% Choose degree')
nx = ChooseDegree(hsv);
end


totbnd=2*tail(nx+1);  % TOTBND = 2*(HSV(NX+1)+HSV(NX+2)+...+HSV(NS))

% 
% Truncate and partition singular vector matrices U,V
u1 = u(1:rh-p,1:nx);
v1 = v(1:ch-m,1:nx);
u2 = u(p+1:rh,1:nx);
u3 = u(rh-p+1:rh,1:nx);
%
sv = sqrt(hsv(1:nx));
invsv = ones(nx,1)./sv;
%

% Kung's formula for the reduced model:
%
% The following is equivalent to UBAR = PINV(U)*[U2; 0]:
ubar = u1'*u2;
a = pinv(u1'*u1)*ubar.*(invsv*sv');   % A = INV(DIAG(SV))*UBAR*DIAG(SV)
b = (sv*ones(1,m)).*v1(1:m,:)';  % B = DIAG(SV)*V1(1:M,:)'
c = u1(1:p,:).*(ones(p,1)*sv');  % C = U1(1:P,:)*DIAG(SV)

%sys = ss(a,b,c,d,-1);
bs=size(b);
cs=size(c);
sys.a=a;
sys.b=b;
sys.c=c;
sys.d=zeros(cs(1),bs(2));


% ----------------------------------------------------------------  myImpulse
function imp = myImpulse(sys,n)
[a,b,c,d] = ssdata(sys);
imp = zeros(n,1);
imp(1) = d;
Ab = b;
for k=2:n
 imp(k) = c*Ab;
 Ab = a * Ab;
end


