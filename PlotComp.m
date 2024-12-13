function PlotComp(S,max_w)

if (nargin<2)
    max_w=50;
end
subplot(2,2,1);
for k=[1:2]
    for j=[1:2]
        subplot(2,2,2*(k-1)+j);
        v=squeeze(S.value(k,j,:));
        plot(real(v),imag(v),'b');
        if (isfield(S,'comp'))
            p=S.comp{k,j};
            w=[linspace(max(S.freq),max_w,1000),linspace(-max_w,min(S.freq),1000)];
            vv=polyval(p,exp(i*L2C(w)));
            hold on;
            plot(real(vv),imag(vv),'r');
            hold off;
        end
    end
end 
% $Id: PlotComp.m,v 1.3 2002/08/30 16:05:00 fseyfert Exp $ 
