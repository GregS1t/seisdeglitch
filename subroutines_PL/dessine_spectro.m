function dessine_spectro(wt,ifig,istep,iaxe,dtw,winwidth,w,wdg,datat)

load noise.mat

%%%%%% plot spectrogram
few=1/dtw;
winwidthw=floor(winwidth/dtw);
overlapw=floor(winwidthw*.5);
nfftw=winwidthw;
hrs=(wt(1)-floor(wt(1)))*24;
if hrs > 24
	hrs = hrs-24;
end
%hrs=(wt(1)-floor(wt(1)))*24-24;

%%%%% spectromgram frequency limit
f1=0.005;
f2=5;


figure(ifig+(istep-1)*2)
subplot(2,3,iaxe)
[S,F1,T,P]=spectrogram(wdg,winwidthw,overlapw,nfftw,few);
T=T+hrs*3600;
F10=F1;
P0=P;
[FF,TT]=meshgrid(F1,T/3600);
pcolor(TT,FF,20*log10(sqrt(P))')
set(gca,'Yscale','log','FontSize',13)
shading interp
shading flat
colormap(gca,'jet')
xlabel(gca,'Time (LMST)')
ylabel(gca,' Hz ')
ylim([f1 f2])
colorbar
if istep == 1 & iaxe == 1
title('First deglitch Z')
elseif istep == 1 & iaxe == 2
title('First deglitch NS')
elseif istep == 1 & iaxe == 3
title('First deglitch EW')
elseif istep == 2 & iaxe == 1
title('Second deglitch without offset')
end
hold on
if iaxe == 1
	loglog(fvel,20*log10(nvel),'k-','LineWidth',2);
else
	loglog(fvel,20*log10(nvel/1.15^2),'k-','LineWidth',2);
end

subplot(2,3,3+iaxe)
[S,F1,T,P]=spectrogram(w,winwidthw,overlapw,nfftw,few);
T=T+hrs*3600;
F1raw=F1;
Praw=P;
[FF,TT]=meshgrid(F1,T/3600);
pcolor(TT,FF,20*log10(sqrt(P))')
set(gca,'Yscale','log','FontSize',13)
shading interp
shading flat
colormap(gca,'jet')
xlabel(gca,'Time (LMST)')
ylabel(gca,' Hz ')
ylim([f1 f2])
colorbar
if istep == 1 & iaxe == 1
title('Raw')
elseif istep == 2 & iaxe == 1
title('Second deglitch with offset')
end
% maximum of scale 
[a1,i1]=min(abs(FF(:,1)-f1));
[a1,i2]=min(abs(FF(:,1)-f2));
cmax=max(max(20*log10(sqrt(P(i1:i2,:)))));
cmax=90.;
subplot(2,3,3+iaxe)
if datat == 'A'
caxis([-200 -140])
else
caxis([1 cmax])
end
subplot(2,3,iaxe)
if datat == 'A'
caxis([-200 -140])
else
caxis([1 cmax])
end

hold on
if iaxe == 1
	loglog(fvel,20*log10(nvel),'k-','LineWidth',2);
else
	loglog(fvel,20*log10(nvel/1.15^2),'k-','LineWidth',2);
end


figure(ifig+iaxe+2)
%%% spectres horaires
Nhr=6;
NNhr=length(P0(1,:));
NHR=floor((TT(end,1)-TT(1,1))*Nhr)+1;
nhr(1:NHR)=0.;
nfhr=length(P0(:,1));
PHR(1:nfhr,1:NHR)=0.;
PHRaw(1:nfhr,1:NHR)=0.;
for i=1:NNhr
ihr=floor((TT(i,1)-TT(1,1))*Nhr)+1;
nhr(ihr)=nhr(ihr)+1;
PHR(1:nfhr,ihr)=PHR(1:nfhr,ihr)+P0(1:nfhr,i);
PHRaw(1:nfhr,ihr)=PHRaw(1:nfhr,ihr)+Praw(1:nfhr,i);
end
for ihr=1:NHR
PHR(1:nfhr,ihr)=sqrt(PHR(1:nfhr,ihr)/nhr(ihr));
PHRaw(1:nfhr,ihr)=sqrt(PHRaw(1:nfhr,ihr)/nhr(ihr));
FHR(1:nfhr,ihr)=F1(1:nfhr);
IHR(1:nfhr,ihr)=(ihr-1)/Nhr+TT(1,1);
Ratio(1:nfhr,ihr)=PHR(1:nfhr,ihr).*power(PHRaw(1:nfhr,ihr),-1);
end
if istep == 1
title('First deglitch')
else
title('Second deglitch')
end
subplot(1,3,1)
plot3(IHR,FHR,PHR,'g*')
hold on
plot3(IHR,FHR,PHRaw,'r*')
set(gca,'Yscale','log','FontSize',13,'Zscale','log')
ylim([5e-3 1.])
subplot(1,3,2)
plot3(IHR,FHR,Ratio,'g*')
set(gca,'Yscale','log','FontSize',13)
subplot(1,3,3)
[a1,i1]=min(abs(FHR(:,1)-.1));
[a1,i2]=min(abs(FHR(:,1)-.5));
semilogy(IHR(i1,:),Ratio(i1,:),'k-','LineWidth',2)
hold on
semilogy(IHR(i2,:),Ratio(i2,:),'r-','LineWidth',2)
xlabel('Time')
ylabel('Amplitude DU/Hz^{1/2}')
legend('10sec','2sec')

figure(ifig+6)
subplot(2,3,iaxe)
loglog(FHR,PHRaw,'.')
hold on
if datat == 'A'
ylim([1e-11 1e-6])
else
ylim([1e-1 1e5])
end
xlim([1e-3 1])
subplot(2,3,iaxe+3)
loglog(FHR,PHR,'.')
hold on
if datat == 'A'
ylim([1e-11 1e-6])
else
ylim([1e-1 1e5])
end
xlim([1e-3 1])

F1=F10;
P=P0;
figure(ifig+7)
subplot(2,3,iaxe)
statistics;
title('Deglitched data')
subplot(2,3,iaxe+3)
F1=F1raw;
P=Praw;
statistics;
title('Raw data')






