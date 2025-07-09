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
function [val,g,v,ecm,ecc,ed,Hes]=EvaluateLagrangian(lmult,L,t_flag)
% Evaluates the langrangian at a dual point lmult. L is a structure
% containing all neccessary parameters

if (nargin<3)
	t_flag=0;
end
lmult=reshape(lmult,length(lmult),1);
P=L.P;
Pc=L.Pc;
r=100/L.norme_f2;
H=r*L.H;
y=L.y;
yc=L.yc;
f=r*L.f;
tP=transpose(conj(P));
tH=transpose(conj(H));
tPP=tP*P;
tHH=tH*H;
switch (L.cflag)
case 1
    % Relaxed problem - minimize only causality
    v=(tHH)\(-tH*f);
    val=norm(H*v+f)^2;
    g=0;
    ed=-1;
case 2,
    % Partial constrained problem - minimize causality under mod. constraint
    lmult_m=lmult;
    if (size(lmult)~=size(yc)) 
        error('Lagrange multiplier should be of size length(yc)');
    end
    tPc=transpose(conj(Pc));
    tPcPc=tPc*diag(lmult_m)*Pc;
    F=inv(tPcPc+tHH);	
    v=F*(-tH*f);
    g=abs(Pc*v).^2-abs(yc).^2;
    G=[tPc*diag(Pc*v)];
    tG=transpose(conj(G));
    Hes=-2*real(tG*F*G);
    ee=norm(H*v+f)^2+transpose(lmult)*(g);
    val=ee;
    ed=-1;
case 3,
    % Full constrained problem - minimize dist. to data under causality constraint and mod. constraint
    sl=length(lmult);
    if (sl~=size(yc)+1) 
        error('Lagrange multiplier should be of size length(yc)+1');
    end
    lmult_m=lmult(1:sl-1);    
    lmult_c=lmult(sl);
    tPc=transpose(conj(Pc));
    tPcPc=tPc*diag(lmult_m)*Pc;
    F=inv(tPP+tPcPc+lmult_c*tHH);
    v=F*(tP*y-lmult_c*tH*f);
    e=P*v-y;
    g=abs(Pc*v).^2-abs(yc).^2;
    g=[g;norm(H*v+f)^2-(L.causalc)^2];
    G=[tPc*diag(Pc*v),tHH*v+tH*f];
    tG=transpose(conj(G));
    Hes=-2*real(tG*F*G);
    ee=transpose(lmult)*(g);
    val=norm(e)^2+ee;
    ed=norm(e);
end;
ecm=max(abs(Pc*v));

%aux=H*v;
%aux(1:10)
%f(1:10)
%theta_aux=linspace(5.6199-2*pi,0.6435,1000);
%y_aux=polyval(transpose(Reverse(v)),1./C2L(theta_aux));
%r*CompFourier(theta_aux,y_aux,[1:10]);

ecc=norm(H*v+f);

%% Hessian test 
if (t_flag)
	kk=length(lmult);
        H_t=[];
	eps_t=0.000000001;
	for j=1:kk
        lmult_t=lmult;
	[val_t,g_t0,v_t,ecm_t,ecc_t,ed_t,Hess_t]=EvaluateLagrangian(lmult_t,L,0);
        lmult_t(j)=lmult_t(j)+eps_t;
	[val_t,g_t1,v_t,ecm_t,ecc_t,ed_t,Hess_t]=EvaluateLagrangian(lmult_t,L,0);
        H_t=[H_t,(g_t1-g_t0)/eps_t];
	end
	fprintf('Error in hessian %f %% Error in sym: %f %%',norm(Hes-H_t)/norm(Hes)*100,norm(Hes-transpose(Hes))/norm(Hes)*100);
end	
        
 

