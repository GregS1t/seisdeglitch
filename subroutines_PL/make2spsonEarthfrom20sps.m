%addpath('../aux_functions');
addpath('/Users/greg/Dropbox (IPGP)/vbb_data_for_philippe/new_tools/aux_functions')

%%%% read data
dirDROP = '/Users/greg/Dropbox (IPGP)/'
%dirDROP='/Users/philippe/Dropbox (IPGP)/';
%dirDROP='/Users/lognonne/Dropbox (IPGP)/';
%%%% read data
%dirU=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/20sps/hg/vbbu'];
%dirV=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/20sps/hg/vbbv'];
%dirW=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/20sps/hg/vbbw'];
dirU=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/20sps/hg/vbbu'];
dirV=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/20sps/hg/vbbv'];
dirW=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/20sps/hg/vbbw'];

%%%  read data
% dirUm=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbue'];
% dirVm=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbve'];
% dirWm=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbwe'];

dirUm=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/2sps/hg/vbbue'];
dirVm=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/2sps/hg/vbbve'];
dirWm=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/2sps/hg/vbbwe'];

%%%  write data
% dirUw=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbue'];
% dirVw=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbve'];
% dirWw=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbwe'];

dirUw=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/2sps/hg/vbbue'];
dirVw=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/2sps/hg/vbbve'];
dirWw=[dirDROP 'DeglitchPkg_PL_GS/DataLucieRolland/2sps/hg/vbbwe'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  3 componenents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Toverlap=5*60;   % this is the overlap of each sol the sol before
sec_earth=86400;
sec_mars=88775.244147;
dtraw=0.05;
dtrawmars=dtraw*sec_earth/sec_mars;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run1:  182 to 236
%  237 NOK
% run2: 238 to 260
% run3: 262 to 267
% 268 to 287 no data
% run 4 288 to 293
% 294,295, 296 NOK
%  run 4 297:303
% run 5: 306 to 350
iweek=0;
for daystart=793:795
first(1:3)=0;
tlast(1:3)=0.;
for icomp=1:3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if icomp == 1
	dirC=dirU;
	dirCs=dirUm;  % synchronisation on the timing of U for V and W also
	dirCm=dirUm;
	dirCw=dirUw;
elseif icomp == 2
	dirC=dirV;
	dirCs=dirUm;  % synchronisation
	dirCm=dirVm;
	dirCw=dirVw;
elseif icomp == 3
	dirC=dirW;
	dirCs=dirUm;  % synchronisation
	dirCm=dirWm;
	dirCw=dirWw;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% this part shall start at sol 168
dayfin=daystart;
%%%%%% data parameters
daylast=182;
daywrite=182;
dayfirst=182;
%%%%% start 
if daystart -1 >= dayfirst
id0=1;
else
id0=0;
end
id1=0;
for iday=daystart:dayfin
clearvars -EXCEPT iday daystart dayfin id0 id1 daylast daywrite dayfirst dirC dirCm dirCw icomp dirU dirV dirW dirUm dirVm dirWm dirUw dirVw dirWw first tlast Toverlap dtraw dtrawmars dirCs iweek
%%%%%% read the day before and the current day
N0=1;
% day before
sday=num2str(iday-1);
iveille=0;
filedata=[dirC,'/sol',sday,'.mat'];
	if exist(filedata,'file') >= 1
	load(filedata);
	NS=length(t);
	utlmst(N0:N0+NS-1)=tlmst;
	ut(N0:N0+NS-1)=t;
	u(N0:N0+NS-1)=ts;
	N0=N0+NS;
	iveille=1;
end
% current day 
sday=num2str(iday);
clear t tlmst ts
filedata=[dirC,'/sol',sday,'.mat'];
load(filedata);
NS=length(t);
utlmst(N0:N0+NS-1)=tlmst;
ut(N0:N0+NS-1)=t;
u(N0:N0+NS-1)=ts;
N0=N0+NS;
%%%%% detect and flag with NaN the interuption
Ninter=1;
Nfirst(1)=1;
for i=2:length(ut)
if (ut(i)-ut(i-1))*24*3600 > 1.25*dtraw
	Ninter=Ninter+1;
	Nfirst(Ninter)=i;
	Duration(Ninter)=(ut(i)-ut(i-1))*24*3600;
else
	Nlast(Ninter)=i;
end
end
if Ninter > 1
'Warning... this segment ( previous and current day) has interuption'
Ninter-1
end
%%%% update time series
k=1;
i=1;
uti(1)=ut(1);
utlmsti(1)=utlmst(1);
N=length(ut);
while (i < N)
dti=ut(i+1)-uti(k);
dtlmsti=utlmst(i+1)-utlmsti(k);
if dti > 1.1*dtraw/24/3600
 	dti=dtraw/24/3600;
 	dtlmsti=dtrawmars/24/3600;
%	dti=dtilast;
%	dtlmstli=dtlmstlast;
	uti(k+1)=uti(k)+dti;
	utlmsti(k+1)=utlmsti(k)+dtlmsti;
	ui(k+1)=NaN;
	
else
	dtilast=dti;
	dtlmstlast=dtlmsti;
	uti(k+1)=uti(k)+dti;
	utlmsti(k+1)=utlmsti(k)+dtlmsti;
	ui(k+1)=u(i+1);
	i=i+1;
end
k=k+1;
end
%%%%% substitute
ut=uti;
utlmst=utlmsti;
u=ui;
%%%%% read PFO filter %%%%%%%%%%%%%%%
PFO5=load('PFO_div5.txt');
gain5=sum(PFO5);
N5=length(PFO5);
%%%%% read PFO filter %%%%%%%%%%%%%%%
PFO2=load('PFO_div2.txt');
gain2=sum(PFO2);
N2=length(PFO2);
%%%%%% synchronise with data of previous sol
Nstart=floor((N2-1)*.25/.5)+1;
Delay=(N5-1)*.05/2/(24*3600)+(N2-1)*.25/2/(24*3600);
if iveille == 1
	sday=num2str(iday-1);
	if iday == daylast
		filedata=[dirCs,'/sol',sday,'-1.mat'];
	else
		filedata=[dirCs,'/sol',sday,'.mat'];
	end
	clear t ts tlmst
	load(filedata);
	utlmst2spslast=tlmst(end);
	ut2spslast=t(end);
	u2spslast=ts(end);
%%%%%%%
if first(icomp) == 0
	first(icomp) = 1;
	tlast(icomp)=ut2spslast;
end
%%%%%%%% synchronize with data on Mars. This define the next data.
% 	[cc i0]=min(abs(ut-ut2spslast-(.5-Nstart*.5*5)/(24*3600)));
  	[cc i0]=min(abs(ut-ut2spslast+10/60/24));
else
	i0=1;
end
%%%%%%%% synchronize with data on Mars
iend=length(ut);
%%% decimate starting with synchronized data   
clear u2spsd ut2spsd utlmst2spsd u4spsd
u4spsd=decimeFPGA(u(i0:iend),5,PFO5,gain5);
u2spsd=decimeFPGA(u4spsd,2,PFO2,gain2);
ut2spsd=downsample(ut(i0:iend),5*2)-Delay;
utlmst2spsd=downsample(utlmst(i0:iend),5*2)-Delay;
if iveille == 1
	[cc i1]=min(abs(ut2spsd-ut2spslast-(.5+Delay-Toverlap)/(24*3600)));
else
	i1=1;
end
%%%%% restart at correct point
u2sps=u2spsd(i1:end);
ut2sps=ut2spsd(i1:end);
utlmst2sps=utlmst2spsd(i1:end);
%%%%%% check with data of sol
figure(10)
subplot(3,1,icomp)
clear tlmst t ts utlmst2spsm ut2spsm u2spsm
if iday <= daylast
sday=num2str(iday);
filedata=[dirCm,'/sol',sday,'.mat'];
load(filedata);
utlmst2spsm=tlmst;
ut2spsm=t;
u2spsm=ts;
plot(utlmst2sps,u2sps,'r-o',utlmst2spsm,u2spsm,'k-*')
else
plot(utlmst2sps,u2sps,'r-o')
end
hold on
%%%% find start and end of sol
i1=1;
for i=1:length(ut2sps)
if utlmst2sps(i)+Delay < iday-Toverlap/24/3600
	i1=i+1;
end
if utlmst2sps(i)+(Delay) < (iday+1)
	i2=i;
end
end
%%%%%%%%%
clear t ts tlmst
t=ut2sps(i1:i2)';
ts=u2sps(i1:i2)';
tlmst=utlmst2sps(i1:i2)';
%%%%% save interuption information
Ninter2sps=1;
Nfirst2sps(1)=1;
i=1;
while (i < length(t))
if ts(i) ~= ts(i)
	Nlast2sps(Ninter2sps)=i-1;
	i=i+1;
	while (i<length(t) & ts(i) ~= ts(i) )
		i=i+1;
	end
	Ninter2sps=Ninter2sps+1;
	Nfirst2sps(Ninter2sps)=i;
end
i=i+1;
end
Nlast2sps(Ninter2sps)=i;
%%%%%%%%
if Ninter2sps > 1
'Warning... this segment has interuption'
Ninter2sps-1
else
'.......... this segment has no interuption'
end
%%%%%%%%%
if iday >= daywrite
sday=num2str(iday);
filedata=[dirCw,'/sol',sday,'.mat'];
save(filedata,'t','ts','tlmst')
if Ninter > 1 || Ninter2sps > 1
filedata=[dirCw,'/intersol',sday,'.mat'];
save(filedata,'Ninter','Nfirst','Nlast','Ninter2sps','Nfirst2sps','Nlast2sps')
end
end
%%%%% cancel NaN   %%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(ts)
if ts(i) ~= ts(i)
	ts(i) = 0.;
end
end
for i=1:length(u)
if u(i) ~= u(i)
	u(i)=0.;
end
end
%%%%% check timing %%%%%%%%%%%%%%%%%%%%%%%%
%%% filter
[B2,A2]=butter(5,.5/2.,'low');
[B10,A10]=butter(5,.5/10.,'low');
figure(20)
subplot(3,1,icomp)
tsf=filtfilt(B2,A2,ts);
uf=filtfilt(B10,A10,u);
plot(datetime(datevec(t)),tsf,'ro')
hold on
plot(datetime(datevec(ut)),uf,'k')
%%%%%% plot drift and data interruption
figure(30)
if icomp == 1
N=length(t);
%delta=-1.6e-6
delta=-1.55e-6
dt0=0.5*(1+delta);
Drift=t(1)+[0:N-1]*dt0/(24*3600);
dDrift=(t-Drift')*24*3600*1e3;
if tlast(icomp) ~= 0
dtjunc=((t(1)-tlast(icomp)+Toverlap/24/3600)*24*3600-.5)*1e3;
else
dtjunc=0.;
end
plot(datetime(datevec(t)),dDrift,'k')
hold on
ylabel('milliseconds')
%ylim([-30 30])
title('Drift')
plot(datetime(datevec(t(1))),dtjunc,'ro')
for i=1:Ninter2sps
plot(datetime(datevec(t(Nfirst2sps(i)))),dDrift(Nfirst2sps(i)),'c<')
plot(datetime(datevec(t(Nlast2sps(i)))),dDrift(Nlast2sps(i)),'c>')
end
%tlast(icomp)=t(end),
end
end
end
if iweek > 7
	iweek=0
%	close all
else
	iweek=iweek+1;
end
end

