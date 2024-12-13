function SaveS(freqs,S,filename)

n=length(freqs);
si=size(S);
assert(n==si(3),'Size of freq should be the same as number of S-mesurements');


fid=fopen(filename,'w');
assert(fid~=-1,'Can not create file');
for k=1:n
    fprintf(fid,'%f ',freqs(k));
    for m=1:2
        for l=1:2
            fprintf(fid,'%f %f ',real(S(m,l,k)),imag(S(m,l,k)));
        end
    end
    fprintf(fid,'\n');
end