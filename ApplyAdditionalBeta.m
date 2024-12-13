function [S,sys,p]=ApplyAdditionalBeta(S,beta,signe_rule,sys,p,val_anti_diag)

[x,k_zero]=min(abs(S.freq));
% Fix the sign of the imaginary part of S12 in zero 
if (IsFilter(S) & imag(exp(i*beta)*(S.value(1,2,k_zero)))*signe_rule<0 )
    beta=beta+pi;
    if (nargin >3) 
        sys.c(1,:)=-1*sys.c(1,:);
        sys.b(:,1)=-1*sys.b(:,1);
        p{1,2}=-p{1,2};
        p{2,1}=-p{2,1};
    end
end
m=exp(i*beta);
[imax,jmax]=GetSDataSize(S);
for k=1:imax
    for j=1:jmax
        if (k~=j)
            S.value(k,j,:)=m*S.value(k,j,:);
            if (isfield(S,'comp'))
                S.comp{k,j}=m*S.comp{k,j};
            end
        end
    end
end
if (nargin>6) 
    for k=1:imax
        for j=1:jmax
            if (k~=j)
                val_anti_diag(k,j)=m*val_anti_diag(k,j);
            end
        end
    end
end
S.b(1,2)=S.b(1,2)+beta;
S.b(2,1)=S.b(2,1)+beta;