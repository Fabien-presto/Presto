function WriteS(S,filename)
%  WriteS(S,filename)
% S: S-structure to write
% filename: destination filename
% Action
% Write the fields freq and value of S to an ascii file. Each line in
% the ascii file is structured as follows:
% S.freq(i) real(S.value(1,1,i)) imag(S(1,1,i)) real(S(1,2,i)) imag(S(1,2,i)) real(S(2,1,i)) imag(S(2,1,i) real(S(2,2i)) imag(S(2,2,i))

n=length(S.freq);
fid = fopen(filename,'w');
if (fid)
	for k=1:n
		fprintf(fid,'%1.9e ',S.freq(k));
		for m=1:2
			for l=1:2
				fprintf(fid,'%1.9e %1.9e ',real(S.value(m,l,k)),imag(S.value(m,l,k)));
			end
		end
		fprintf(fid,'\n');
	end
	fclose(fid);
else
	error('Cant write on file');
end
    



