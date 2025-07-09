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
%line 96 <example.nw>
function [nD , G] = dataRotate(data, X, allpass)
nD = data;
c = exp(i*X(1));
nD.value(2,1,:) = c * nD.value(2,1,:) ;
nD.value(1,2,:) = c * nD.value(1,2,:) ;

if isfield(nD, 'valinf')
  nD.valinf(2,1,:) = c * nD.valinf(2,1,:) ;
  nD.valinf(1,2,:) = c * nD.valinf(1,2,:) ;
end


if nargout==1, return, end

[n,m,p]   = SSdim(allpass);
[P,M,LL]  = size(nD.value);
[a,b,c,d] = SSdata(allpass);

switch data.type
case 'sample'
 %%%% ****** Warning **************
 %%%% Here it is explicetally supposed that the parematrisation of the all
 %%%% pass factor is done in C and A  - B being computed by projection 
 %%%% There should be a switch exactly equal to the one used in
 %%%% arl2UserFunction
 
 z = exp(i*nD.theta); 
 % Old non vectorised code
 %for k=1:LL, zA{k} = inv(z(k)*eye(size(a))-a); end
 
 %Vectorised code
 aux_sys=ss(a,eye(size(a)),c,zeros(p,n));
 CzA=freqresp(aux_sys,z);
 
 
 WB = zeros(n,M); dWB = zeros(n,M); 
 WD = zeros(p,M); dWD = zeros(p,M); 
 W = zeros(n,n); L = zeros(n,p);
 % Old non vectorised code 
%  for k = 1:LL
%   Fk = nD.value(:,:,k) ; CzA = c*zA{k}; 
%   WB = WB + CzA' * Fk; WD = WD + Fk;
%   dWB = dWB + CzA' * (Fk .* [0 i;i 0]);
%   dWD = dWD + (Fk .* [0 i;i 0]);
%   W  = W  + CzA' * CzA ;
%   L  = L  + CzA';
%  end
 % Vectorised code
 F=nD.value(:,:,:) ;
 WB= sum(MatrixProd3D(MatrixStar3D(CzA),F),3);WD=sum(F,3);
 dWB= sum(MatrixProd3D(MatrixStar3D(CzA),F.*repmat([0 i;i 0],[1,1,LL])),3);
 dWD= sum(F.*repmat([0 i;i 0],[1,1,LL]),3);
 W= sum(MatrixProd3D(MatrixStar3D(CzA),CzA),3);
 L= sum(MatrixStar3D(CzA),3);
 
 
 
 W = W/LL; WB = WB/LL; dWB = dWB/LL; L = L/LL; WD = WD/LL; dWD = dWD/LL;
 X = [W L ; L' eye(p,p)] \ [WB ; WD];
 B = X(1:n,:); D = X(n+1:end,:);
 G =   -2* real(trace(D'*dWD + B'*dWB));

case 'coef'
 % Left allpass factor, best linear estimate of B and D matrices
 B = zeros(n,m); dB = zeros(n,m);
 for k = LL:-1:1
  x = nD.value(:,:,k);
  B = a'* B + c'* x;
  dB = a'* dB + c'* (x .* [0 i; i 0]);
 end

 G = -2*real(trace(B'*dB));
end

