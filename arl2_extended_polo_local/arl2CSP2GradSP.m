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
%line 200 <allpass.nw>
function dy = arl2CSP2GradSP(CSP,DJ,Slist,Ulist,Vlist)

[p, N] = size(CSP.y);
%line 205 <allpass.nw>
Lambda = DJ;
dy = [];

% Loop on degree
for n = N:-1:1
 % Current Schur parameter
 v    = CSP.y(:,n);
 w    = CSP.chart.w(n);
 u    = CSP.chart.u(:,n);
 v2 = v'*v;
 w2 = w'*w;
 Iw2=sqrt(1-w2);
 Iv2=sqrt(1-v2);
 Iwv2=sqrt(1-w2*v2);
 uu=u*u';
 Id=eye(size(uu));
 Zuu=zeros(size(uu));
 Zu=zeros(size(u));

 UU   = Ulist{n};
 VV   = Vlist{n};
 SS   = Slist{n};

 DV = Lambda*UU'*SS';
 DU = SS'*VV'*Lambda;
 II = 1:p+1;
 DV = DV(II,II);
 DU = DU(II,II);

% HERE IS THE TRICKY CODE
 % Only alpha for U
 gv1 =   w2/(Iwv2^3);
 gv2 = - Iw2^2 / ( Iv2 * Iwv2^3);

 dU1 = [Iw2*u , Zuu ; 0 , Iw2*u'];
 dU2 = [ Zu , -w*uu ; w' , Zu'  ];

 alphaU = scal(DU, dU1')*gv1 + scal(DU, dU2')*gv2;

 % alpha and beta for V
 dV1 = [ Iw2*v ,   Zuu  ; 0  , -Iw2*v'];
 if v2 > 0
  dV2 = [ Zu    , v*v'/(v2)  ; 1  ,   Zu'  ];
 else
  dV2 = [ Zu    , Zuu        ; 1  ,   Zu'  ];
 end

 alphaV = scal(DV, dV1)*gv1 + scal(DV, dV2)*gv2;

 beta  = (Iw2/Iwv2)*(DV(1:p,1) - DV(end, 2:end)');

 % Upper right term
 DVUR = DV(1:p,2:end);

 if v2 > 0
  alphaV = alphaV + 2*(1-(Iv2/Iwv2))*scal(DVUR, v*v' )/(v2^2) ;
  beta = beta     - (1-(Iv2/Iwv2))*(DVUR+DVUR')*v/v2;
 end
 dyn = (alphaU + alphaV)*v + beta;

 dy = [dyn , dy];
%line 268 <allpass.nw>
 % Apply (FUV) transform
 Lambda = VV' * Lambda * UU' ;
 Lambda = Lambda(2:end, 2:end);
%line 273 <allpass.nw>
% End loop on degree
end
%line 277 <allpass.nw>
function y = scal(A, B)
y = real(sum(sum(A .* conj(B))));
