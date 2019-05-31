% script to write a wavelet from my 
% convoluted alpha function network activity 
if(idx_1==1)Outp=b1;end %%Icells
if(idx_1==0)Outp=b2;end %%Ecells
% script to write a wavelet from my 
% convoluted alpha function network activity

f_low=0.01;  % lowest frequency
f_high=70;  % highest frequency
f_step=0.1; % frequency step
%%figure
[W,t,fq]=Wavelet_1ch(Outp,fs,f_low,f_high,f_step);
imagesc(t,fq,flipud(abs(W)));
colormap(jet(64))
axis xy;
set(gca,'FontSize',15,'XTick',[10,11,12,13,14,15,16],'YTick',[5,10,15,20,25,30],'ticklength',3.5*get(gca,'ticklength'),'TickDirMode','manual','TickDir','out')
cb=colorbar();
set(cb,'FontSize',10);
xlim([10 16]);
ylim([5 30]);
xlabel('time (s)')
ylabel('Frequency (Hz)')
caxis([0 70]) % change colorbar limits
box off;
in_v=input1;