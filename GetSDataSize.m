function [imax,jmax,n]=GetSDataSize(S)

assert(isfield(S,'value'),'No data in S structure');
[imax,jmax,n]=size(S.value);