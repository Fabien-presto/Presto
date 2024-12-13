function [ss]=CayleyReal(s)

Id=eye(size(s.d));
aux=inv(s.d+Id);
ss.a = s.a - s.b*aux*s.c;
ss.b=i*sqrt(2)*s.b*aux;
ss.c=i*sqrt(2)*aux*s.c;
ss.d=(Id-s.d)*aux;
