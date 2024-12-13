function [M,Mred,R1,R2,Rya]=S2M(S,flag)
% function [M,Mred,R1,R2]=S2M(S,flag)
% Computes the exended coupling matrix as well as the reduced coupling
% matrix of the rational response contained in S and previousely computed
% by presto
% Optional flag={0,1,2} specifies which arrow is wanted or if the transversal for is wanted - either ML1=0 or
% ML1=MSN - flag=2 corresponds to the transversal form


assert(isfield(S,'p') & isfield(S,'q'),'S contains no rational approximant');
if nargin<2 
    flag=0;
end
R=Presto_Guilbert(S.p,S.q);
% Compute realization of admittance response
Ry=CayleyReal(R);
% Computing arrow form
if (flag==2)
    [P,aux]=eig(Ry.a);
    P=P*sqrt(pinv(P.'*P));
    Ry.a = P.'*Ry.a*P;
    Ry.b = P.'*Ry.b;
    Ry.c = Ry.c*P;
    Rya=Ry;
else
    if (length(Ry.a)>1)
        [P,Rya]=ToArrow(Ry,flag);
    else
        Rya=Ry;
    end
    % Adjust sign of loads
    SP=eye(length(Rya.a),length(Rya.a));
    SP(1,1)=sign(real(Rya.b(1,1)));
    SP(end,end)=sign(real(Rya.b(end,2)));
    Rya.a = SP*Rya.a*SP;
    Rya.b = SP*Rya.b;
    Rya.c = Rya.c*SP;
end

M=zeros(length(Rya.a)+2,length(Rya.a)+2);
M(2:end-1,2:end-1)=-Rya.a;
M(1,1)=Rya.d(1,1);
M(end,end)=Rya.d(2,2);
M(end,1)=Rya.d(2,1);
M(1,end)=Rya.d(1,2);
M(2:end-1,1)=Rya.b(1:end,1)*i;
M(2:end-1,end)=Rya.b(1:end,2)*i;
M(1,2:end-1)=Rya.c(1,1:end)*i;
M(end,2:end-1)=Rya.c(2,1:end)*i;
Mred=M(2:end-1,2:end-1);
R1=imag(M(1,2))^2;
R2=imag(M(end,end-1))^2;

