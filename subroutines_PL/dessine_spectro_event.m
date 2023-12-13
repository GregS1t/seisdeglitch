function dessine_spectro_event(wt,iaxe,dtw,winwidth,w,wdg,T0,ah1)


%%%%%% plot spectrogram
few=1/dtw;
winwidthw=floor(winwidth/dtw);
overlapw=floor(winwidthw*.5);
nfftw=winwidthw;
mins=(wt(1)-T0)*24*60;

%%%%% spectromgram frequency limit
f1=0.04;
f2=1;


ah2=subplot(3,3,3+iaxe)
[S,F1,T,P]=spectrogram(wdg,winwidthw,overlapw,nfftw,few);
T=T+mins*60;
F10=F1;
P0=P;
[FF,TT]=meshgrid(F1,T/60);
pcolor(TT,FF,20*log10(sqrt(P))')
set(gca,'Yscale','log','FontSize',13)
shading interp
shading flat
colormap(gca,'jet')
xlabel(gca,'Time (min)')
ylabel(gca,' Hz ')
ylim([f1 f2])
colorbar
title('Deglitched data')
caxis([-20 60])
xlim([-10 60])

ah3=subplot(3,3,6+iaxe)
[S,F1,T,P]=spectrogram(w,winwidthw,overlapw,nfftw,few);
T=T+mins*60;
F1raw=F1;
Praw=P;
[FF,TT]=meshgrid(F1,T/60);
pcolor(TT,FF,20*log10(sqrt(P))')
set(gca,'Yscale','log','FontSize',13)
shading interp
shading flat
colormap(gca,'jet')
xlabel(gca,'Time (min)')
ylabel(gca,' Hz ')
ylim([f1 f2])
colorbar
title('Glitch')
caxis([-20 60])
xlim([-10 60])

pos3 = get(ah3,'Position');
pos2 = get(ah2,'Position');
pos1 = get(ah1,'Position');
pos2(3) = pos1(3);
pos3(3) = pos1(3);
set(ah2,'Position',pos2)
set(ah3,'Position',pos3)






