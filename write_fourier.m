rs='TEST_';
fid=fopen([rs,'11'],'w');
fprintf(fid,'%g %g \n',[real(Sre2.fourier{1,1}(1:end));imag(Sre2.fourier{1,1}(1:end))]);
fclose(fid);
rs='TEST_';
fid=fopen([rs,'12'],'w');
fprintf(fid,'%g %g \n',[real(Sre2.fourier{1,2}(1:end));imag(Sre2.fourier{1,2}(1:end))]);
fclose(fid);
rs='TEST_';
fid=fopen([rs,'21'],'w');
fprintf(fid,'%g %g \n',[real(Sre2.fourier{2,1}(1:end));imag(Sre2.fourier{2,1}(1:end))]);
fclose(fid);
fid=fopen([rs,'22'],'w');
fprintf(fid,'%g %g \n',[real(Sre2.fourier{2,2}(1:end));imag(Sre2.fourier{2,2}(1:end))]);
fclose(fid);
    
    
     
% $Id: write_fourier.m,v 1.3 2010/10/12 15:25:42 fseyfert Exp $ 
