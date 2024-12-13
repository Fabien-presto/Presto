function PlotErrors(S)

w=squeeze(S.freq);
assert(isfield(S,'p')& isfield(S,'q'),'No rational approximation attached to the S-structure');
sd=FreqResp(S.p,S.q,w);
subplot(2,2,1);
for k=1:2
    for j=1:2
        subplot(2,2,2*(k-1)+j);
        plot(S.freq,abs(squeeze(sd(k,j,:))-squeeze(S.value(k,j,:))));
        %plot(S.freq,abs(sd(k,j)));
        %plot(S.freq,abs(squeeze(S.value(k,j,:))));
    end
end
