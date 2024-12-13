function [H]=CreateBlocHankel(blocs)

n=length(blocs);
assert( n>0 & mod(n,1)==0,'No blocs or no odd number of blocs');
[l,m]=size(blocs{1});
nn=(n+1)/2;
H=zeros(l*nn,m*nn);
for k=1:nn
    for j=1:k
        assert(prod([l,m]==size(blocs{k})),'Blocs should be of same size'); 
        H((k-j)*l+1:(k-j+1)*l,(j-1)*m+1:(j)*m)=blocs{k};        
    end
end
offset=l*nn;
for k=nn+1:n
    kk=k-nn;
    for j=1:nn-kk
        assert(prod([l,m]==size(blocs{k})),'Blocs should be of same size'); 
        H(offset+(kk-j-1)*l+1:offset+(kk-j)*l,j*m+1:(j+1)*m)=blocs{k};        
    end
end