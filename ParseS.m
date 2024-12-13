function [S]=ParseS(filename,format_flag)
% [S]=ParseS(filename,format_flag)
% filename: string specifiying then name of the file to parse
% Action
% Parses an ascii file into an S structure.
% If format_flag (default) is ri then each ligne is suppposed to be in
% following format
% format flag=ri
% freq(i) real(S11(i)) imag(S11(i)) real(S12(i)) imag(S12(i)) real(S21(i)) imag(S21(i) real(S22(i)) imag(S22(i)
% format_flag=mad
% freq(i) abs(S11(i)) angle(S11(i)) abs(S12(i)) angle(S12(i)) abs(S21(i))angle(S21(i) abs(S22(i)) angle(S22(i)
% where angle is expressd in degrees
% format_flag=mar 
% Same as mad but angle in radiant 
% format_flag=dbd 
% Magnitude in db and angle in degree 

if nargin==1 
    format_flag='ri';
end

switch format_flag
    case 'ri'
        [freq,re11,im11,re12,im12,re21,im21,re22,im22]=textread(filename,'%f %f %f %f %f %f %f %f %f ','commentstyle','shell');
        
        S.freq=freq;
        S.value=zeros(2,2,length(S.freq));
        S.value(1,1,:)=re11+i*im11;
        S.value(1,2,:)=re12+i*im12;
        S.value(2,1,:)=re21+i*im21;
        S.value(2,2,:)=re22+i*im22;
        
    case 'mad'
        [freq,m11,a11,m12,a12,m21,a21,m22,a22]=textread(filename,'%f %f %f %f %f %f %f %f %f ','commentstyle','shell');
        S.freq=freq;
        S.value=zeros(2,2,length(S.freq));
        r=pi/180;
        S.value(1,1,:)=m11.*exp(i*a11*r);
        S.value(1,2,:)=m12.*exp(i*a12*r);
        S.value(2,1,:)=m21.*exp(i*a21*r);
        S.value(2,2,:)=m22.*exp(i*a22*r);
     case 'mar'
        [freq,m11,a11,m12,a12,m21,a21,m22,a22]=textread(filename,'%f %f %f %f %f %f %f %f %f ','commentstyle','shell');
        S.freq=freq;
        S.value=zeros(2,2,length(S.freq));
        r=1;
        S.value(1,1,:)=m11.*exp(i*a11*r);
        S.value(1,2,:)=m12.*exp(i*a12*r);
        S.value(2,1,:)=m21.*exp(i*a21*r);
        S.value(2,2,:)=m22.*exp(i*a22*r); 
     case 'dbd'
        [freq,m11,a11,m12,a12,m21,a21,m22,a22]=textread(filename,'%f %f %f %f %f %f %f %f %f ','commentstyle','shell');
        S.freq=freq;
        S.value=zeros(2,2,length(S.freq));
        r=pi/180;
        S.value(1,1,:)=dbtomag(m11).*exp(i*a11*r);
        S.value(1,2,:)=dbtomag(m12).*exp(i*a12*r);
        S.value(2,1,:)=dbtomag(m21).*exp(i*a21*r);
        S.value(2,2,:)=dbtomag(m22).*exp(i*a22*r);
end
end
function [mag]=dbtomag(db)
mag=10.^(db/20);
end
% $Id: ParseS.m,v 1.5 2002/09/03 13:24:33 fseyfert Exp $ 
