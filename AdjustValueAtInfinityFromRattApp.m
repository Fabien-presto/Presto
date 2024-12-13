function S=AdjustValueAtInfinityFromRattApp(S)

[imax,jmax]=GetSDataSize(S);
if (~isfield(S,'p') | ~isfield(S,'q'))
    error('No rational approximation present');
end
a=[];
for k=1:min(imax,jmax)
    if length(S.q)==length(S.p{k,k})
        a(k)=-angle(S.p{k,k}(1)/S.q(1));
    end
end
S=ApplyDiagonalRotation(S,a);
        
        