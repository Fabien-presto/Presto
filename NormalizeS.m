function S=NormalizeS(S,fmin,fmax);
% S=NormalizeS(S,fmin,fmax)
% Normalize the frequency of S by sending the passband [fmin,fmax] to [-1,1] via a
% linear transformation 
fmin=min(fmin,fmax);
fmax=max(fmin,fmax);
a=2/(fmax-fmin);
b=-a*(fmax+fmin)/2;
S.freq=a*S.freq+b;