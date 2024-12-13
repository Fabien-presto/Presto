function [b]=IsFilter(S)

[im,jm]=GetSDataSize(S);
b=(im==2 & jm==2);