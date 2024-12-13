function [S,sys]=AllSteps(S,n,zeros_at_inf,solver_flag,plot_flag)
% [Sres,sys]=AllSteps(S,n,zeros_at_inf,solver_flag,plot_flag)
% S: S measurements
% n: target MacMillan degree
% zeros_at_inf (optional): annulation order at infty for S12 and S21
% solver_flag (optional) : rational approx. engine used - {'RARL2','HARL2'}
% plot_flag (optional) : {0,1}
% Action
% Proceeds to the following steps:
%          - CompensateDelayAndFreqShift (and build completion)
%          - CompFourierCoeffs
%          - RatApp (launches rational approximation engine)
% The result is S enriched with the completion, the compensations, the rational approximation
% The output sys is a bootstrap strucutre than can be used in further
% identifications as starting point - typically in the function RatAppLocal
 

assert(nargin>1,'No McMillan degree given');
S=ReshapeS(S);
% Get Default values and Constants
DVC=GetDefaultValuesAndConstants(S);

if (nargin<5)
    plot_flag=DVC.AS.plot_flag;
end
if (nargin<4 | isempty(solver_flag))
    solver_flag=DVC.AS.solver_flag;
end
if (nargin<3 | isempty(zeros_at_inf))
    zeros_at_inf=DVC.AS.zeros_at_inf;
end
if (nargin==3 & ischar(zeros_at_inf))
    solver_flag=DVC.AS.solver_flag;
    zeros_at_inf=DVC.AS.zeros_at_inf;
end
if (plot_flag)
    clf;
    PlotS(S);
    RaiseAndDraw;
end
fprintf('Compensate delay and freq. shift - compute completion \n');
tic;
S=CompensateDelayAndFreqShift(S,zeros_at_inf);
fprintf('Time to compensate delays and shift:%f\n',toc);
if (plot_flag)
    PlotS(S);
    RaiseAndDraw;
    %pause;
end
fprintf('Compute  fourier coeff. - take stable part\n');
tic;
S=CompFourierCoeffs(S);
fprintf('Time to compute fourier coeffs:%f\n',toc);
if (plot_flag)
    PlotS(S);
    RaiseAndDraw;
end
fprintf('Compute  rational approximation \n');
t_rat_app=tic();
[S,sys]=RatApp(S,n,DVC.RA.divide_by_z_minus_one,solver_flag);
S=AdjustValueAtInfinityFromRattApp(S);
fprintf('Time to compute rat. app. %f\n',toc(t_rat_app));
if (plot_flag)
    PlotS(S);
    RaiseAndDraw;
end
 
% $Id: AllSteps.m,v 1.9 2010/10/12 15:25:42 fseyfert Exp $ 
