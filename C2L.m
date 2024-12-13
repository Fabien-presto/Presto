function [w]=C2L(theta)
assert(norm(imag(theta))<=1e-6,'Angles. are not real');
w=imag(moeb(exp(i*theta))); 
% $Id: C2L.m,v 1.2 2002/08/30 16:05:00 fseyfert Exp $ 
