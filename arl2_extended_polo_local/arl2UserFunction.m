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
%line 6 <userF.nw>
function [err, sys, gradJ] = arl2UserFunction(allpass , data)

switch nargout
 case 1, needGrad = 0;
 case 2, needGrad = 0;
 case 3, needGrad = 1;
 case 4, needGrad = 1;
 otherwise, error('Bad output args count')
end

[a,b,c,d] = SSdata(allpass);
[n,m,p]   = SSdim(allpass);
arl2Assert(p==m);

switch data.type
case 'sys'
   
%line 51 <userF.nw>
TheSys      = data.value;
[N, M, P]   = SSdim(TheSys);
if min(M, P) ~= p ,   error('arl2UserFunction: dimensions mismatch'), end
[A,B,C,D] = SSdata(TheSys);
i1 = 1:N; i2 = N+1:N+n;
Ag = [A zeros(N,n) ; zeros(n,N) a];  Dg = zeros(size(D));
%line 59 <userF.nw>
if P < M 
 % Left allpass factor, best linear estimate of B matrix
 Cg = [C, c]; WC = dlyap(Ag' , Cg'*Cg);
 W21 = WC(i2,i1);  W22 = WC(i2,i2);
 arl2Assert(arl2IsSmall(norm(W22 - eye(size(W22)))) ,  'Not allpass ?' );
 x =  W21*B; Bg = [B ; -x]; 
 err = SSnH2(SScreate(Ag,Bg,Cg,Dg));
 sys = SScreate(a,x,c,D);
 % grad(J) wrt A and C
 if needGrad
   WB = dlyap(Ag , Bg*Bg'); W12 = WB(i1,i2);
   DJA = W21*A*W12; DJC =  C*W12;
   gradJ = 2*[ zeros(p+n, p), [DJC; DJA] ];
 end
%line 75 <userF.nw>
else
 Bg = [B; b]; WB = dlyap(Ag , Bg*Bg');
 % Right allpass factor, best linear estimate of C matrix, 
 W12 = WB(i1,i2); W22 = WB(i2,i2);
 arl2Assert(arl2IsSmall(norm(W22 - eye(size(W22)))) ,  'Not allpass ?' );
 x =  C*W12; Cg = [C , -x]; 
 err = SSnH2(SScreate(Ag,Bg,Cg,Dg));
 sys = SScreate(a,b,x,D);
 % grad(J) wrt A and B
 if needGrad
   WC = dlyap(Ag' , Cg'*Cg);  W21 = WC(i2,i1);
   DJA = W21*A*W12; DJB =   W21*B;
   gradJ = 2*[zeros(p, p+n) ; [DJB, DJA] ];
 end
end
%line 23 <userF.nw>
case 'sample'
   
%line 93 <userF.nw>
TheFunct = data.value;
[P,M,LL] = size(TheFunct);
if min(M, P) ~= p,   error('Dimensions mismatch'), end
z = exp(i*data.theta); 



%line 100 <userF.nw>
%%% No testing on the dimension for the moment - because incompatible with
%%% dataRotate function where parametrization is always supposed to be in A
%%% and C 
if true
%if P <= M
    %for k=1:LL, zA{k} = inv(z(k)*eye(size(a))-a); end
    %aux_sys=ss(a,eye(size(a)),c,zeros(p,n));
    %CzA=freqresp(aux_sys,z);    
    if (needGrad)
        aux_sys=ss(a,eye(size(a)),eye(size(a)),zeros(size(a)));
        zA=freqresp(aux_sys,z);
        cc=repmat(c,[1,1,LL]);
        CzA=MatrixProd3D(cc,zA);
    else
        aux_sys=ss(a,eye(size(a)),c,zeros(p,n));
        CzA=freqresp(aux_sys,z);
    end
    % Left allpass factor, best linear estimate of B matrix
    WB = zeros(n,M); W = zeros(n,n); L = zeros(n,p); WD = zeros(p,M);
    F0 = zeros(M,M);
    % Old code without vecotrisation
    %  for k = 1:LL
    %   Fk = TheFunct(:,:,k); %CzA = c*zA{k};
    %   WB = WB + CzA(:,:,k)' * Fk; WD = WD + Fk;
    %   W  = W  + CzA(:,:,k)' * CzA(:,:,k) ;
    %   L  = L  + CzA(:,:,k)';
    %   F0 = F0 + Fk'*Fk;
    %  end
    
    % Vectorised code
    WB=sum(MatrixProd3D(MatrixStar3D(CzA),TheFunct),3);WD=sum(TheFunct,3);
    W=sum(MatrixProd3D(MatrixStar3D(CzA),CzA),3);
    L=sum(MatrixStar3D(CzA),3);
    F0=sum(MatrixProd3D(MatrixStar3D(TheFunct),TheFunct),3);
    
    W = W/LL; WB = WB/LL; L = L/LL;  WD =  WD/LL ; F0 = F0/LL;
    X = [W L ; L' eye(p,p)] \ [WB ; WD];
    B = X(1:n,:); D = X(n+1:end,:);
    sys = SScreate(a,B,c,D);
    err = real(trace(F0 - D'*WD - B'*WB));
    
    if needGrad
        %aux_sys=ss(a,B,eye(size(a)),zeros(n,m));
        %zAB=freqresp(aux_sys,z);
        BB=repmat(B,[1,1,LL]);
        zAB=MatrixProd3D(zA,BB);
        %aux_sys=ss(a,B,c,D);
        %CzABD=freqresp(aux_sys,z);
        DD=repmat(D,[1,1,LL]);
        CzABD=MatrixProd3D(cc,zAB)+DD;
        
        DJC = zeros(p,n);
        DJA = zeros(n,n);
        % Old non vectorised code
        %   for k = 1:LL
        %    %zAB = zA{k}*B; CzA = c*zA{k};
        %    ek = TheFunct(:,:,k) - CzABD(:,:,k);
        %    errSignal(:,:,k) = 2*ek/LL;
        %    DJC = DJC - 2*ek*zAB(:,:,k)';
        %    DJA = DJA - 2*CzA(:,:,k)'*ek*zAB(:,:,k)';
        %   end
        % Vectorised code
        e=TheFunct-CzABD;
        errSignal=2*e/LL;
        DJC=-2*sum(MatrixProd3D(e,MatrixStar3D(zAB)),3);
        DJA=-2*sum(MatrixProd3D(MatrixProd3D(MatrixStar3D(CzA),e),MatrixStar3D(zAB)),3);
        
        
        
        
        
        gradJ = [ zeros(p+n, p), [DJC; DJA] ]/LL;
    end
%line 130 <userF.nw>
else
    if (needGrad)
        aux_sys=ss(a,eye(size(a)),eye(size(a)),zeros(size(a)));
        zA=freqresp(aux_sys,z);
        bb=repmat(b,[1,1,LL]);
        zAB=MatrixProd3D(zA,bb);
    else
        aux_sys=ss(a,b,eye(size(a)),zeros(n,m));
        zAB=freqresp(aux_sys,z);
    end
    
    % Right allpass factor, best linear estimate of C matrix
    CW = zeros(P,n); W = zeros(n,n); L = zeros(n,m); DW = zeros(P,m);
    F0 = zeros(P,P);
    %  Old vectorised code
    %  for k = 1:LL
    %   Fk = TheFunct(:,:,k); zAB = zA{k}*b;
    %   CW = CW + Fk  * zAB'; DW = DW +  Fk;
    %   W  = W  + zAB * zAB';
    %   L  = L  + zAB;
    %   F0 = F0 + Fk*Fk';
    %  end
    
    % Vectorised code
    CW=sum(MatrixProd3D(TheFunct,MatrixStar3D(zAB)),3);DW=sum(TheFunct,3);
    W=sum(MatrixProd3D(zAB,MatrixStar3D(zAB)),3);
    L=sum(zAB,3);
    F0=sum(MatrixProd3D(TheFunct,MatrixStar3D(TheFunct)),3);
    
    W = W/LL; CW = CW/LL; L = L/LL; DW = DW/LL; F0 = F0/LL;
    X =  [CW  , data.D] / [W L ; L' eye(m,m)] ;
    C = X(:,1:n); D = X(:,n+1:end);
    sys = SScreate(a,b,C,D);
    err = real(trace(F0 - D*DW' - C*CW'));
    
    if needGrad
        CC=repmat(C,[1,1,LL]);
        CzA=MatrixProd3D(CC,zA);
        DD=repmat(D,[1,1,LL]);
        CzABD=MatrixProd3D(CzA,bb)+DD;
        DJB = zeros(n,m);
        DJA = zeros(n,n);
        % Old non-vectorised code
        %   for k = 1:LL
        %    zAB = zA{k}*b; CzA = C*zA{k}; ek = TheFunct(:,:,k) - CzA*b - D;
        %    errSignal(:,:,k) = 2*ek/LL;
        %    DJB = DJB - 2*CzA'*ek;
        %    DJA = DJA - 2*CzA'*ek*zAB';
        %   end
        
        % Vectorised Code
        e=TheFunct-CzABD;
        errSignal=2*e/LL;
        DJB=-2*sum(MatrixProd3D(MatrixStar3D(CzA),e),3);
        DJA=-2*sum(MatrixProd3D(MatrixProd3D(MatrixStar3D(CzA),e),MatrixStar3D(zAB)),3);
        
        gradJ =[zeros(p, p+n) ; [DJB, DJA] ]/LL;
    end
    %line 160 <userF.nw>
end

%line 25 <userF.nw>
case 'coef'
   
%line 165 <userF.nw>
TheCoeff = data.value;
[P,M,LL] = size(TheCoeff);
if min(M, P) ~= p,   error('Dimensions mismatch'), end
%line 170 <userF.nw>
%%% No testing on the dimension for the moment - because incompatible with
%%% dataRotate function where parametrization is always supposed to be in A
%%% and C 
if true
%if P <= M 
 % Left allpass factor, best linear estimate of B and D matrices
 B = zeros(n,M);  F0 = zeros(M,M); 

 for k = LL:-1:1
  beta(:,:,k) = B;
  x =TheCoeff(:,:,k);
  B = a'* B + c'* x;
  F0 = F0 +   x' * x;
 end
 D   = data.valinf;
 sys = SScreate(a,B,c,D);
 err = real(trace(F0 - B'*B));

 if needGrad
  DJC = zeros(p,n); 
  DJA = zeros(n,n); 
  for k = LL:-1:1
   DJC = DJC * a' - 2*TheCoeff(:,:,k)*B';
   DJA = DJA * a' - 2*beta(:,:,k)*B';
  end

  gradJ = [ zeros(p+n, p), [DJC; DJA] ];
 end
%line 196 <userF.nw>
else
 % Right allpass factor, best linear estimate of C matrix
 C = zeros(P,n);  F0 = zeros(P,P); 
 for k = LL:-1:1
  gamma(:,:,k) = C;
  x =TheCoeff(:,:,k);
  C = C * a' +  x * b' ;
  F0 = F0 +   x * x';
 end
 D   = data.valinf;
 sys = SScreate(a,b,C,D);
 err = real(trace(F0 - C*C'));

 if needGrad
  DJB = zeros(n,m); 
  DJA = zeros(n,n); 
  for k = LL:-1:1
   DJB = a' * DJB  - 2*C'*TheCoeff(:,:,k);
   DJA = a' * DJA  - 2*C'*gamma(:,:,k);
  end

  gradJ =[zeros(p, p+n) ; [DJB, DJA] ];
 end
%line 221 <userF.nw>
end

 

%line 27 <userF.nw>
otherwise
   error('Bad data type');
end
