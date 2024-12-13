function [theta]=L2C(w)
assert(norm(imag(w))<=1e-6,'Freq. are not real');
theta=mod(imag(log(moeb(i*w))),2*pi); 
% $Id: L2C.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
