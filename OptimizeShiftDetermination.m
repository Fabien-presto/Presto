function [beta]=OptimizeShiftDetermination(fouriers,n)

m=100;
beta_range=linspace(-pi/50,pi/50,m);

crits=[];
mfouriers=fouriers;
for r=beta_range
    edelta=exp(i*r);
    for k=1:2*50+1
        fouriers{1,2}(k)=mfouriers{1,2}(k)*edelta;
        fouriers{2,1}(k)=mfouriers{2,1}(k)*edelta;
    end
    for k=1:2*50+1
        for l=1:2
            for m=1:2
                f(l,m)=fouriers{l,m}(k+1);
            end
        end
        blocs{k}=f;
    end
    H=CreateBlocHankel(blocs);
    sl=svds(H,n+1)
    crits=[crits,sl(n+1)];
end
[critmin,imin]=min(crits);
beta=beta_range(imin);
plot(beta_range,crits);
        