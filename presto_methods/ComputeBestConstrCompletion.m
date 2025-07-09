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
function [v]=ComputeBestConstrCompletion(w,S,w_lim,module_c,causal_c,n,n_cp,n_cf,iso_flag)
% Solves several hierarchically contrained problems for the completion problem
% Optimisation criterion: l2 distance to data
% Always feasible contraint: modulus uper bound of the modulus of
% completion
% Other contraint: Distance to causal function of completed data

eps=0.01;
n_spline_points=5000; 

if (module_c < 0)
	error('Modulus contraint is <0');
end

theta=L2C(w);
theta_lim=L2C(w_lim);
theta1=2*pi-theta_lim;
theta2=theta_lim;
ind=find(theta<theta1 | theta>theta2);
if (length(ind)<n+3) 
    warning('Very few data points to perform approximation at that order');
end
x=w(ind);
y=S(ind);
theta_max=min(theta)+2*pi;
theta_min=max(theta);
thetab=linspace(min(theta),max(theta),n_spline_points);
% Divide S by z-1 if using L^2 isometry 
if (iso_flag)
	z=transpose(exp(i*theta));
	S=S./(z-1);
end
yb=spline(theta,S,thetab);
f=(FourierCoeffs(thetab,yb,[1:n_cf]));
% Build Structure for Lagrangian computation
sf.f=f;
sf.theta1=theta_min;
sf.theta2=theta_max;
if (iso_flag)
	% S is already divided by z-1 but there is still a factor
	% to match the norm on the plan
	sf.norm_f2=sqrt(SimpleQuad(theta,abs(S).^2)/(pi));
else
	sf.norm_f2=sqrt(SimpleQuad(theta,abs(S).^2)/(2*pi));
end
% Build constraint vector
xc=ComputeTchebyPoints(1/w(1),1/w(end),n_cp);
yc=spline([1/w(1),1/w(end),0],[abs(S(1)),abs(S(end)),module_c],xc);

Ls=ComputeStructureForLagrangian(1./x,y,xc,yc,sf,causal_c,n,1,iso_flag);

% Change of basis in the lagrangian - choose Cheby poly for the completion
[Ls,Pc]=ReconditionLagrangian(Ls,1/w(1),1/w(end));

l=[];
Ls.cflag=1;
[val,g,v,ecm,ecc]=EvaluateLagrangian(l,Ls);
% Try one shot without contraints
l=zeros(size(xc)+[0,1]);
Ls.cflag=3;
[val,g,v,ecm,ecc]=EvaluateLagrangian(l,Ls);
if (ecm<module_c+eps & ecc<causal_c+eps )
    fprintf('No constraint optimisation needed (constr. completion)\n');	
    v=Reverse(transpose(Pc*v));	
    return;
else
    % With constraint on modulus
    %options = optimset('GradObj','on','Hessian','on','HessUpdate','bfgs','DerivativeCheck','off','Diagnostics','off','LargeScale','on','Display','off','TolFun',1e-8,'TolX',1e-8);
    options = optimset('GradObj','on','Algorithm','trust-region-reflective','Hessian','on','DerivativeCheck','off','Diagnostics','off','LargeScale','on','Display','off','TolFun',1e-8,'TolX',1e-8);
    Ls.cflag=2;
    l=ones(size(xc));
    fprintf('Start partial constraint optimisation \n');
    for num_try=1:3
        [l]=fmincon(@OptLagrangian,l,[],[],[],[],zeros(size(l)),[],[],options,Ls);
        [val,g,v,ecm,ecc]=EvaluateLagrangian(l,Ls);
        if (sum(g>eps))
               options=optimset('TolX',optimget(options,'TolX')/5);
               fprintf('Improve numeric precision \n');
        else
               break;
        end
    end   
    if (ecc>causal_c-eps)
        fprintf('Causal constraint unfeasable (>%1.3f) - take best possible \n',ecc);
	    v=Reverse(transpose(Pc*v));
        return;
    else
        % With all constraints
        options = optimset(options,'TolX',1e-8);
        fprintf('Start full constraint optimisation \n');
        Ls.cflag=3;
        l=[l,0];
        [val,g,v,ecm,ecc,ed]=EvaluateLagrangian(l,Ls);
        for num_try=1:3
            [l,fval,exitflag,output,lambda,grad,hessian]=fmincon(@OptLagrangian,l,[],[],[],[],zeros(size(l)),[],[],options,Ls);
            [val,g,v,ecm,ecc,ed]=EvaluateLagrangian(l,Ls);
            if (sum(g>eps))
                options=optimset(options,'TolX',optimget(options,'TolX')/5);
                fprintf('Improve numeric precision \n');
            else
                break;
            end
        end
        v=Reverse(transpose(Pc*v));
    end
end
 

