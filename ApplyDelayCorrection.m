function [Sout]=ApplyDelayCorrection(Sin,Sref)
% [Sout]=ApplyDelayCorrection(Sin,Sref)
% Applies the phase correction Sref.a and Sref.b to the data Sin

a=Sref.a;
b=Sref.b;
Sout=Sin;
for k=1:2
    for j=1:2
        Sout.value(k,j,:)=squeeze(Sout.value(k,j,:)).*(exp(a(k,j)*1i*Sout.freq));
        Sout.value(k,j,:)=Sout.value(k,j,:)*exp(1i*b(k,j)); 
    end
end
if (~isfield(Sin,'a'))
    Sout.a=zeros(2,2);
end
if (~isfield(Sin,'b'))
    Sout.b=zeros(2,2);
end
Sout.a=Sout.a+a;
Sout.b=Sout.b+b;

 
% $Id: ApplyDelayCorrection.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
