function [LOCM,NLOC,PEAK,datauf,datavf,datawf]=locate_glitch(dataut,datau,datavt,datav,datawt,dataw,dt)
disp("In locate_glitch function")
[B,A]=butter(2,[.001*dt/.5 ,.4*dt/.5]);         %%%% filter remove low freq
%[B,A]=butter(2,[.0001 ,.4]);         %%%% filter remove low freq
datauf=(filtfilt(B,A,datau));            %%%
datavf=(filtfilt(B,A,datav));             %%%
datawf=(filtfilt(B,A,dataw));           %%%
dataufa=abs(datauf);
datavfa=abs(datavf);
datawfa=abs(datawf);

% First plot
figure(41)
sgtitle("Pick search")
subplot(1,2,1)
plot(datetime(datevec(dataut)), dataufa, 'r', "DisplayName", "u");
hold on
plot(datetime(datevec(datavt)), datavfa, 'b', "DisplayName", "v");
plot(datetime(datevec(datawt)), datawfa, 'g', "DisplayName", "w");
legend
grid on

kk=0;
% QUESTION
% Le même vecteur de localisation de glitch est utilisé pour les 3 axes  ????
for i=2:length(datau)-1
    if  (( datauf(i)>datauf(i-1) & datauf(i)>datauf(i+1) )  || ( datauf(i)<datauf(i-1) & datauf(i)<datauf(i+1) ))
        kk=kk+1;                            %%%% peak index
        peak(kk)=abs(datauf(i));
        loc(kk)=i;                             %%%% location
        loct(kk)=dataut(i);		         %%%% time
    end
end

for i=2:length(datav)-1 
     if  ( ( datavf(i)>datavf(i-1) & datavf(i)>datavf(i+1) )  || ( datavf(i)<datavf(i-1) & datavf(i)<datavf(i+1) ))
       kk=kk+1;                                     	%%%% peak index
       peak(kk)=abs(datavf(i));
       loc(kk)=i;                             %%%% location
	loct(kk)=datavt(i);		       %%%% time
    end
end
for i=2:length(dataw)-1
     if  ( ( datawf(i)>datawf(i-1) & datawf(i)>datawf(i+1) )  || ( datawf(i)<datawf(i-1) & datawf(i)<datawf(i+1) ))
       kk=kk+1;                                     	%%%% peak index
       peak(kk)=abs(datawf(i));
       loc(kk)=i;                             %%%% location
	loct(kk)=datawt(i);		       %%%% time
    end
end
%%

plot(datetime(datevec(loct)), peak, 'r*');

[peakm, iloc]=sort(peak);
locm=loc(iloc);
loctm=loct(iloc);
%threshold=.8;
%n1=floor((1-threshold)*kk);
n1=1;

%LOCM=locm(n1:kk);
%NLOC=length(LOCM);
%PEAK=peakm(n1:kk);

% Bypass the sort of the vectors by amplitude
PEAK  = peak(n1:kk);
LOCM = loc(n1:kk);
NLOC=length(LOCM);

%%%%%%
subplot(1,2,2)
plot(datetime(datevec(dataut)),datauf,'k', "DisplayName", "u");
hold on
plot(datetime(datevec(dataut(locm(n1:kk)))),datau(locm(n1:kk)),'rO', "DisplayName", "pick U");
plot(datetime(datevec(datavt)),datavf,'k',  "DisplayName", "v");

plot(datetime(datevec(datavt(locm(n1:kk)))),datav(locm(n1:kk)),'bO', "DisplayName", "pick V");
plot(datetime(datevec(datawt)),datawf,'k',  "DisplayName", "u");

plot(datetime(datevec(datawt(locm(n1:kk)))),dataw(locm(n1:kk)),'gO', "DisplayName", "pick W");
title("Peaks identified for each axis")
legend
grid on
end
