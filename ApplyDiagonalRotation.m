function S=ApplyDiagonalRotation(S,a)

[imax,jmax]=GetSDataSize(S);

for k=1:imax
    for j=1:jmax
        if (k==j)
            S.value(k,k,:)=exp(i*a(k)).*S.value(k,k,:);
            if (isfield(S,'p'))
                S.p{k,k}=exp(i*a(k)).*S.p{k,k};
            end
            S.b(k,k)=angle(exp(i*S.b(k,k))*exp(i*a(k)));
        else
            if (k<=jmax)
                S.value(k,j,:)=exp(i*a(k)/2).*S.value(k,j,:);
                if (isfield(S,'p'))
                    S.p{k,j}=exp(i*a(k)/2).*S.p{k,j};
                end
                S.b(k,j)=angle(exp(i*S.b(k,j))*exp(i*a(k)/2));
            end
            if (j<=imax)
                S.value(k,j,:)=exp(i*a(j)/2).*S.value(k,j,:);
                if (isfield(S,'p'))
                    S.p{k,j}=exp(i*a(j)/2).*S.p{k,j};
                end
                S.b(k,j)=angle(exp(i*S.b(k,j))*exp(i*a(j)/2));
            end
        end
    end
end

    