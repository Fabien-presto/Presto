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
function [S]=CompensateDelayAndFreqShift(S,zeros_at_inf,plotflag)
% [Sres]=CompensateDelayAndFreqShift(S,deg_ratio,zeros_at_inf,plotflag)
% S: S measurements
% deg_ratio: degree of the expected rational model
% zeros_at_inf (optional): annulation order at infty for S12 and S21
% plot_flag (optional): {0,1}
% Action
% For S11 ans S22 tries to fit measurments at infty with an 1/s serie so as to be able
% to compensate delays.
% For S12 and S21 determines the freq shift which leeds to the less dissipativ system (in a sens).
% Important: See GetDefaultValuesAndConstants.m for Constants used by the algorithm
% The result is S enriched with the completion and the various compensations

Test_date=0;
Limit_date='19-May-2005';

if (Test_date==1)
    date_now=datenum(clock);
    delta_date=date_now-datenum(Limit_date);
    if (delta_date > 0 )
        display(['Sorry, temporary license has expired on ', Limit_date]);
        return;
    else
        display(['Temporary license expires on ',Limit_date]);
    end
end

S=ReshapeS(S);
DVC=GetDefaultValuesAndConstants(S);
[imax,jmax]=GetSDataSize(S);

% Magic numbers for delay determination
iso_flag=DVC.CDAFS.iso_flag;
delay_range=DVC.CDAFS.delay_range;
if (length(DVC.CDAFS.poly_order)==1)
    for k=1:imax
        for j=1:jmax
           poly_order(k,j)=DVC.CDAFS.poly_order;
        end
    end
else
    poly_order=DVC.CDAFS.poly_order;
end
omega_lim=DVC.CDAFS.omega_lim;
error_lim=DVC.CDAFS.error_lim;
sign_in_zero_trans_12=DVC.CDAFS.sign_in_zero_trans_12;
strategy=DVC.CDAFS.strategy;
n_cf=DVC.CDAFS.number_of_fourier_coeffs;
% Hack to be harmonised
modulus_factor=DVC.CDAFS.modulus_factor(1,1);
causal_bound=DVC.CDAFS.causal_bound(1,1);
if (strcmp(strategy,'border_match'))
    fprintf('Strategy chosen for delay det. is: border_match\n');
else
    fprintf('Strategy chosen for delay. det is: causal \n');
end
if (nargin<3)
    zeros_at_inf=DVC.CDAFS.zeros_at_inf;
end
if (nargin<4)
        plotflag=DVC.CDAFS.plot_flag;
end
% Create a cell array for polynomial expansion
S.comp=cell(imax,jmax);
S.a=zeros(imax,jmax);
S.b=zeros(imax,jmax);
for (k=1:imax)
    for (j=1:jmax)
          Skj=transpose(squeeze(S.value(k,j,:)));
          if ((k==j))
              if (strcmp(strategy,'border_match'))
                    [S.comp{k,j},S.a(k,j),errmin]=DetDelay(S.freq,Skj,delay_range,poly_order(k,j),0,omega_lim);
	                if(errmin>error_lim)
                    fprintf('Warning - Error is sup. to limit in delay and comp. det. \n');
                    end
              else
                  [S.a(k,j),svdmin]=ComputeBestCausalCompensation(S.freq,Skj,delay_range,poly_order(k,j),n_cf,omega_lim,modulus_factor);
              end
              fprintf('Delay [%d,%d]:%f \n',k,j,S.a(k,j));
              S.value(k,j,:)=Skj.*exp(i.*S.freq*S.a(k,j));
              Skj=transpose(squeeze(S.value(k,j,:)));
              % Start completion determination for (1,1) and (2,2)
              module_c=modulus_factor*max(abs(Skj(:)));
              causal_c= causal_bound;
              n_cp= DVC.CDAFS.number_of_control_points;
              [S.comp{k,j}]=ComputeBestConstrCompletion(S.freq(:),Skj,omega_lim,module_c,causal_c,poly_order(k,j),n_cp,n_cf,iso_flag);
              S.b(k,j)=S.comp{k,j}(length(S.comp{k,j}));
              S.b(k,j)=imag(log(abs(S.b(k,j))/S.b(k,j)));
              for l=1:jmax
                  if (k~=l)                      
                    S.b(k,l)=S.b(k,l)+S.b(k,j)/2;
                  end
              end
              for l=1:imax
                  if (j~=l)
                      S.b(l,j)=S.b(l,j)+S.b(k,j)/2;
                  end
              end
              S.value(k,j,:)=exp(i*S.b(k,j)).*S.value(k,j,:);
              S.comp{k,j}=exp(i*S.b(k,j)).*S.comp{k,j};
          end
    end
end
% Apply b compensation to diagonal terms - applying this b preservs innerness 
for k=1:imax
    for j=1:jmax
        if (k~=j)
            S.value(k,j,:)=exp(i*S.b(k,j)).*S.value(k,j,:);
        end
    end
end

% Compensate and compute completion of diagonal terms
for (k=1:imax)
    for (j=1:jmax)
        if (k~=j)
            a=0;
            if k<=jmax
                a=a+S.a(k,k)/2;
            end
            if (j<=imax)
                a=a+S.a(j,j)/2;
            end
            Skj=transpose(squeeze(S.value(k,j,:)));
            Skj=Skj.*exp(i*a*S.freq);
            S.value(k,j,:)=Skj;
            [S.comp{k,j},S.a(k,j),errmin]=DetDelay(S.freq,Skj,[0],poly_order(k,j)+zeros_at_inf,zeros_at_inf,omega_lim,0);
            S.a(k,j)=a;
            %fprintf('[%d,%d]: Err. 1/s serie app. (order %d) %1.3f%% \n',k,j,poly_order,errmin*100);
        end
    end
end

% Determine anti-diagonal b's (called Beta)
if (IsFilter(S))
    [beta]=DetBeta(S);
    S=ApplyAdditionalBeta(S,beta,sign_in_zero_trans_12);
end
if (plotflag)
    PlotS(S);
end

% Store isometric information
S.comp_iso_flag=iso_flag;
    
       
            
 

