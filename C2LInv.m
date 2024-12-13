function [w]=C2LInv(theta)
assert(norm(imag(theta))<=1e-6,'Angles. are not real');
w=imag(-1*moeb(exp(i*theta),1)); 