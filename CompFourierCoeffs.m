function [S]=CompFourierCoeffs(S)
% [Sres]=CompFourierCoeffs(S)
% S: S measurements
% Action
% Computes Fourier coeff. of S on the circle (using L^{infty} isometrie)
% S should contain information about completion at infty
% See GetDefaultValuesAndConstants.m for Constants used by the algorithm (for ex. number of coeff.)


S=ReshapeS(S);
[imax,jmax]=GetSDataSize(S);
DVC=GetDefaultValuesAndConstants(S);
% Magic numbers for fourrier coeffs computation
n_comp=DVC.CFC.n_comp;
n_spline=DVC.CFC.n_spline;
n_fourier=DVC.CFC.n_fourier;
assert(isfield(S,'comp') &isfield(S,'comp_iso_flag') ,'S has no completion to infinity');
iso_flag=S.comp_iso_flag;
if (isfield(S,'fourier'))
    rmfield(S,'fourier');
end
if (isfield(S,'ratio'))
    rmfield(S,'ratio');
end
S.fourier=cell(imax,jmax);
S.ratio=zeros(imax,jmax);
for j=1:imax
    for k=1:jmax
        %S.fourier{j,k}=CompScalarFourierCoeffs(S.freq,squeeze(S.value(j,k,:)),S.comp{j,k},[0:n_fourier],n_comp,n_spline);
        %c_neg=CompScalarFourierCoeffs(S.freq,squeeze(S.value(j,k,:)),S.comp{j,k},[-n_fourier,-1],n_comp,n_spline);
        S.fourier{j,k}=CompScalarFourierCoeffs(S.freq,squeeze(S.value(j,k,:)),S.comp{j,k},[-1*n_fourier:0],n_comp,n_spline);
        S.fourier{j,k}(:)=S.fourier{j,k}(n_fourier+1:-1:1);
        c_pos=CompScalarFourierCoeffs(S.freq,squeeze(S.value(j,k,:)),S.comp{j,k},[1:n_fourier],n_comp,n_spline);
        S.fourier_pos{j,k}=c_pos;    
	if (iso_flag)	
		r=norm(DivideByZminusOneScalar(c_pos))/sqrt(SimpleQuad(S.freq,abs(squeeze(S.value(j,k,:))).^2)/(2*pi));
        else
		r=norm(c_pos)/sqrt(SimpleQuad(L2C(S.freq),abs(squeeze(S.value(j,k,:))).^2)/(2*pi));
	end
	fprintf('S(%d,%d) ratio unstable/norme_data = %f%%\n',j,k,r*100);
        S.ratio(j,k)=r;
    end
end
 
% $Id: CompFourierCoeffs.m,v 1.6 2002/09/09 15:47:14 fseyfert Exp $ 
