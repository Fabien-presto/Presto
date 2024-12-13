function [P, SS] = ToArrow(S,flag)
% [P, SS] = ToArrow(S,flag): Puts the system S in general arrow form - flag={0,1} drives the
% determination of MSN anf ML1 couplings - for flag=0 ML1 is set to zero
% while for flag=0 MSN=ML1 is imposed - 0  is default
if nargin<2 
    flag=0;
end
% Test des dimensions
B=S.b;
A=S.a;
[Blig,Bcol]=size(B);
if  Bcol ~= 2  
 mesg=sprintf(...
   'B (%d x %d) : bad sizes - must be \n',...
   Blig, Bcol);
 error(mesg);
end
% Test de l'orthogonalit�
% If transpose(B)*B is diagonal, then alpha=0 V=unit(B)
[v1, v2, alpha] = orthog(B,flag);

% Construction de la r�alisation
n = Blig;
%if (imag(transpose(v1)*B(:,1))<0)
%    v1=-v1;
%end;
%if (imag(transpose(v2)*B(:,2))<0)
%    v2=-v2;
%end;
P = [ v1 zeros(n,n-2) v2 ] ; 
x = v1;
for j=2:n-1
 x = A * x;
 y = x;
 x = unit(x - P * (transpose(P)*x) );
 if (imag(transpose(x)*y)<0)
     x = -x;
 end;
 P(:,j) = x;
end
% Construction de la r�alisation
SS.a = transpose(P)*A*P;
SS.b = transpose(P)*B;
SS.c = S.c*P;
SS.d = S.d;

function  [v1, v2, alpha] = orthog(B,flag)
 v1 = unit(B(:,1));
 v2 = unit(B(:,2));
 alpha = scal(v1,v2);
 if (flag)
    x = alpha/( 1 + sqrt(1-alpha^2) );
    vv = v1;
    v1 = unit(v1 - x*v2);
    v2 = unit(v2 - x*vv);
 else
     v1=unit(v1-alpha*v2);
 end
function y = unit(x)
 y = x / sqrt(sum(x.*x));
function p=scal(u,v)
 p=sum(u.*v); 
