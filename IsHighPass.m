function [b]=IsHighPass(S)
b=sign(S.freq(1)*S.freq(end))>=0;