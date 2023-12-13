function [unew,uglitch,udrift,udiff,udg,u0,...
    vnew,vglitch,vdrift,vdiff,vdg,v0,...
    wnew,wglitch,wdrift,wdiff,wdg,w0,...
    uuznew,uuzdrift,NGLITCH,Var0,Var1,val,val0,...
    idx,Amp,Delay,Green,ND,im,g0,ND2,NG1,NG2,...
    dgmin,dgmax,NZEROALL,spredictor,jdx,AmpSpike,reducA, Var]...
    =deglitch5_full_3axis(ut,utlmst,udg,u,...
    vt,vtlmst,vdg,v,wt,wtlmst,wdg,w,uuzdg,...
    uuz,recomp,epsi,NS,dt,Uz,Up,Vz,Vp,Wz,Wp,...
    NZERO,valthreshold,dtthreshold,tfthreshold,...
    dgmin,dgmax,istep,iday,icorrect,ideriv,igg,iequal);

%%%%%%%%
trend=u(1)+(ut-ut(1))/(ut(end)-ut(1))*(u(end)-u(1));
udg0=udg;
u0=u;
%%%%%%%%
trend=v(1)+(vt-vt(1))/(vt(end)-vt(1))*(v(end)-v(1));
vdg0=vdg;
v0=v;
%%%%%%%%
trend=w(1)+(wt-wt(1))/(wt(end)-wt(1))*(w(end)-w(1));
wdg0=wdg;
w0=w;
%%%%%%%%
trend=uuz(1)+(ut-ut(1))/(ut(end)-ut(1))*(uuz(end)-uuz(1));
uuzdg0=uuzdg;
uuz0=uuz;

%%%%%%%%%
N=length(w0);
if ideriv == 1
    Nderiv=2;
else
    Nderiv=3;
end

%%%% take the w axis
clear LOCM NLOC
%[LOCM,NLOC,PEAK]=locate_glitchLMST(wt,wtlmst,w0,dt,.01,.1,6,17);
%[LOCM,NLOC,PEAK]=locate_glitch(ut,u0,dt);
%[LOCM,NLOC,PEAK]=locate_glitch(ut,uuz0,dt);
[LOCM,NLOC,PEAK]=locate_glitch_3axis(ut,u0,vt,v0,wt,w0,dt);

% Save the value to chech order of the peaks
%save output_locate_glitch_3_axis.mat PEAK  NLOC  LOCM

[unew,uglitch,udrift,udiff,uhp,ustack,ustack0,...
    vnew,vglitch,vdrift,vdiff,vhp,vstack,vstack0,...
    wnew,wglitch,wdrift,wdiff,whp,wstack,wstack0,...
    uuznew,uuzdrift,NGLITCH,Var0,Var1,val,val0,idx,Amp,Delay,...
    Green,ND,im,g0,ND2,NDS,nstack,NG1,NG2,NZEROALL,spredictor,jdx,AmpSpike,reducA, Var]...
    =deglitch7_3axis(ut,utlmst,udg,u,vt,vtlmst,vdg,v,wt,wtlmst,wdg,w,...
    uuzdg,uuz,recomp,epsi,NS,dt,Uz,Up,Vz,Vp,Uz,Up,NZERO,...
    valthreshold,dtthreshold,tfthreshold,LOCM,NLOC,Nderiv,istep,iday,icorrect,igg,iequal);

%%%%%
udg=u0-unew;            % remove the glitch from the detrended
vdg=v0-vnew;            % remove the glitch from the detrended
wdg=w0-wnew;            % remove the glitch from the detrended
uuzdg=uuz0-uuznew;            % remove the glitch from the detrended

for iaxe=1:3

figure(300)
if iaxe == 1
	gcolor='k';
	s=u;
	sdg=udg;
	st=ut;
	shp=uhp;
	sglitch=uglitch;
	sdrift=udrift;
elseif iaxe == 2
	gcolor='b';
	s=v;
	sdg=vdg;
	st=vt;
	shp=vhp;
	sglitch=vglitch;
	sdrift=vdrift;
elseif iaxe == 3
	gcolor='g';
	s=w;
	sdg=wdg;
	st=wt;
	shp=whp;
	sglitch=wglitch;
	sdrift=wdrift;
end

subplot(3,1,1)
type=' 5 mHz High pass';
plot(datetime(datevec(st)),shp,gcolor,'LineWidth',1)
hold on
plot(datetime(datevec(st)),sglitch,'r','LineWidth',1)
plot(datetime(datevec(st)),sdrift,'r--','LineWidth',1)

title([ type 'and glitches'])

%%%%%%%
subplot(3,1,2)

title('Raw and deglitched')
plot(datetime(datevec(st)),s,'-r','LineWidth',1)
hold on
plot(datetime(datevec(st)),sdg, gcolor,'LineWidth',1)
dgmax=max([max(sdg),dgmax]);
dgmin=min([min(sdg),dgmin]);
dgmin=dgmin-0.2*abs(dgmin);
dgmax=dgmax+0.2*abs(dgmax);

subplot(3,1,3)

tau=Delay(:,iaxe);
plot(datetime(datevec(st(idx(1:NGLITCH(iaxe),iaxe)))),tau(1:NGLITCH(iaxe)),'*k')
hold on
%xlim([datetime(datevec(st(1))),datetime(datevec(st(end)))]);
title('Spike delay')

figure(1300)
subplot(4,1,4)
plot(datetime(datevec(st)),s,'r')
hold on
if istep > 1
plot(datetime(datevec(st)),sdg,gcolor,'LineWidth',2)
end
title('Raw and deglitches')

end
end
