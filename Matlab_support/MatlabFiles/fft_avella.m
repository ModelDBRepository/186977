%file to comput the fourier spectrum of my signal 
%using the pwelch function.
Fs=fs;
win_PSD=5*Fs;
noverlap_PSD=[];
nFFT=2.^10;
[spd,f]=pwelch(Outp,win_PSD,noverlap_PSD,nFFT,Fs); %spd=spectral power
                                                  %arbitrary units/Hz

%%figure                                                  
axis([0 60 0 75])
plot(f,spd.^0.5);   %we plot here the fourier amplitude
%plot(f,spd.^0.5);   %to not over enphasyse the peaks.
set(gca,'FontSize',13);
xlabel('Frequency (Hz)')
ylabel('Spectral amplitude (Arbitrary Units/Hz)')                                    
ylim([0 20])%%ylim([0 40])
xlim([0 25])
%%grid