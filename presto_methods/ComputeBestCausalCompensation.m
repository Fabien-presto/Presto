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
function [delay,emin,v]=ComputeBestCausalCompensation(w,y,d_range,poly_order,n_cf,w_lim,modulus_factor)
% Computes best delay compensation by looking for best causal completion 
% delay is the optimal delay found and ecmin is a measure of the non causal part

causal_bound=3; % Bound about causality in % which is required to be a local minimum according to continuity considerations

n_spline_points=max(4000,n_cf+1);
number_control_points=3*(poly_order+1);

w=reshape(w,length(w),1);
y=reshape(y,length(y),1);

theta=L2C(w);
max_thetas=max(theta);
min_thetas=min(theta);
% Going to the equaly spaced circle because of the use of the fft
delta=2*pi/n_spline_points;
min_thetas=ceil(min_thetas/delta)*delta;
max_thetas=floor(max_thetas/delta)*delta;
thetas=[min_thetas:delta:max_thetas+delta/2];
ys=spline(theta,y,thetas);
ws=C2L(thetas);

theta_max=min_thetas+2*pi;
theta_min=max_thetas;
thetas_complementary=[theta_min+delta:delta:theta_max-delta/2];
H=GenerateMoebMatrix(theta_min,theta_max,n_cf,poly_order+1);
tH=transpose(conj(H));
%tHH=tH*H;
%invtHH=inv(tHH);
ecs=[];


L.B=H;
L.BB=(L.B)'*L.B;

max_mod=modulus_factor*max(abs(ys(1)),abs(ys(length(ys))));

w_cp=ComputeTchebyPoints(1/w(1),1/w(end),number_control_points);
L.D=GeneratePowerMatrix(w_cp,poly_order);
L.mc=transpose((spline([1/ws(1),1/ws(end),0],[abs(ys(1)),abs(ys(end)),max_mod],w_cp)));

%Change basis for comp
[L,Pc]=ReconditionCausalLagrangian(L,1/w(1),1/w(end));


% Initialisation of lagrange coeffs 
% 0 is best when the degree of the completion is low but for higher degrees
% a strictly positive value
% is needed (otherwise the completion problem is ill posed, matrix singular

if (poly_order<8)
    l=zeros(number_control_points,1);
else
    l=ones(number_control_points,1);
end

% Finds a lower bound for causality contraints -
L.cflag=1;
norm_y2=sqrt(SimpleQuad(theta,abs(y).^2)/(2*pi));

k_counter=0;
d_range=sort(d_range);
d_range=d_range(end:-1:1);
ecc=[];
i_min=1;
%xplot=linspace(1/ws(1),1/ws(end),100);
%clf;
%hold on;
%plot(w_cp,L.mc,'*r');

for r=d_range
    yc=ys.*exp(i*r*ws);
    % Compute associated fourier coeffs
    if (mod(k_counter,10)==0) 
        fprintf('*');
    end
    k_counter=k_counter+1;
    L.f=transpose(MyLocalFFT(min_thetas,n_spline_points,yc,[1:n_cf]));
    [v,ectest,l,mini_error]=MyLocalOptimiser(l,L);
    v=Pc*v;
    %yplot=polyval(v(end:-1:1),xplot);
    %plot(xplot,abs(yplot));
    yp=[polyval(v(end:-1:1),1/ws(1)),polyval(v(end:-1:1),1/ws(end))];
    e1=abs(yp(1)-yc(1));
    e2=abs(yp(end)-yc(end));
    ecs=[ecs,max([e1,e2,mini_error])];
    ecc=[ecc,100*ectest/norm_y2];
    % Finds first (true !)local minimum starting from biggest delays
    if (k_counter>1)
        if (ecs(k_counter)>1.2*ecs(i_min) & ecc(k_counter)<causal_bound )
            break;
        end
        if (ecs(k_counter)<ecs(i_min))
            i_min=k_counter;
            v_min=v;
        end
    end
end
%yplot=polyval(v_min(end:-1:1),xplot);
%plot(xplot,abs(yplot),'g');
%pause;
%clf;
%plot(d_range(1:k_counter),100*ecs,'r',d_range(1:k_counter),ecc,'b');
%plot(d_range(1:k_counter),ecs);
%pause;
%ecc



%[emin,i_min]=min(ecs);
delay=d_range(i_min);
emin=ecs(i_min);

%r=d_range(i_min);
%yc=ys.*exp(i*r*ws);
%L.f=transpose(MyLocalFFT(min_thetas,n_spline_points,yc,[1:n_cf]));
%[v,ectest,l]=MyLocalOptimiser(l,L);
%v=v(end:-1:1)
%wp=linspace(1/ws(1),1/ws(end),1000);
%clf;
%plotc(yc,'.',polyval(v,wp));
%pause;


function [f]=MyLocalFFT(theta,n,y,range)

m=length(y);
yp=[y(2:m),y(1)];
ym=[y(m),y(1:m-1)];
yt=y/2+(yp+ym)/4;
yt(1)=(y(1)+y(2))/4;
yt(m)=(y(m)+y(m-1))/4;
f=fft(yt,n)/(n);
f=f.*exp(-theta*i*[0:n-1]);
f=f(range+1);

function [f]=MyLocalFFT_bis(theta,y,range)

m=length(y);
yp=[y(2:m),y(1)];
ym=[y(m),y(1:m-1)];
yt=y/2+(yp+ym)/4;
f=fft(yt)/(m);
f=f.*exp(-theta*i*[0:m-1]);
f=f(range+1);

function [v,err,l,mini_error]=MyLocalOptimiser(l,L)

eps=0.001;
lb=zeros(size(l));

%options = optimset('GradObj','on','Hessian','on','HessUpdate','bfgs','DerivativeCheck','off','Diagnostics','off','LargeScale','on','Display','off','TolFun',1e-8,'TolX',1e-8);   
options = optimset('GradObj','on','Algorithm','trust-region-reflective','Hessian','on','DerivativeCheck','off','Diagnostics','off','LargeScale','on','Display','off','TolFun',1e-8,'TolX',1e-8);
for num_try=1:3
    [l]=fmincon(@OptCausalLagrangian,l,[],[],[],[],lb,[],[],options,L);
    [val,v,g,Hvoid,err]=EvaluateCausalLagrangian(l,L,0);
    mini_error=sqrt(sum(abs(g)))/length(g);
    if (mini_error>eps)
        options=optimset(options,'TolX',optimget(options,'TolX')/5);
        %fprintf('+');
    else
        break;
    end
end

