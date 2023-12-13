%%%%% version du 20 Juillet 2020
%%%%  version pour evènement: iextract = 1  pour sol, iextract = 0
%%%% read data

%% DEPENDANCE : 
% DATA:
%    données VBB à 2Hz : vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg
%    réponse instrum   : ./metadata/* 
% FUNCTIONS : 
%   read_resp_v2
%   deglitch5_full_3axis
%       - locate_glitch_3axis
%       - deglitch7_3axis
%           - makeglitch4
%               - decimeFPGA
%           - shift
%   FourierleastsquarewithNaN
%   dessine_spectro
%   equalize
%   remove_equalizing
%   remove_tfvitdt2
%   extractanderase0

clear all 
%%
%% Read input the data 

% Source of the data shouldn't be changed

% Path to Dropbox
%dirDROP = '/Users/deraucourt/Dropbox (IPGP)/'
 dirDROP = '/Users/greg/Dropbox (IPGP)/';
 
% Path to 2 sps VEL High Gain VBB data for the 3 components, U, V, W
dirU=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbue'];
dirV=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbve'];
dirW=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbwe'];

%Dir to the FIR filter files
dir_FIR                 = [dirDROP  'DeglitchPkg_PL_GS/FIR/'];
PFO5=load([dir_FIR 'PFO_div5.txt']);      % FIR coefs to decimate by 5 
PFO2=load([dir_FIR 'PFO_div2.txt']);      % FIR coefs to decimate by 2

% Dir to the metadata
dir_resp               = [dirDROP  'DeglitchPkg_PL_GS/metadata/'];
resp.vbb_vel_hg   = [dir_resp  'RESP.XB.ELYSE.02.BHU_new'];
resp.vbb_pos_hg  = [dir_resp  'RESP.XB.ELYSE.00.LMU'];
resp.sp_hg           = [dir_resp  'RESP.XB.ELYSE.65.EHU'];


%% Primary set of parameters
 
%   iseed 0/1 ( sortie seed aussi)   
%   irot 0/1 deglitch en U,V,W ou X,Y,Z  
%   ifilter = 0 ( doit rester a zero)   
%   ideriv = 0 ( reste a zero)
%   Nhar :  Nbr d'harmique du detrend   
%   istep : nombre de de deglitchage ( 1 ou 2)  
%   epsi = 0 ( peut etre par la suite different de zero)
%   iequal=1 ( 0 pas bien) egalize les TF a celle de U
%   icorrect=0 ( 1 possible )
%

iseed=0;
icorrect=0;
ifilter=0;
irot=1;
ideriv=0;
Nhar=12;
istep=2; % number of steps of deglitching
iextract=0;
valthreshold=0.02;
%valthreshold=0.0;
epsi=0.;
iequal=1;

solstart = 82;
solend  = 82;

%% Output directories with deglichted files in both MSEED and MATformat
% WARNING : Must be customized with 
% DG is the glitches selection for threshold = 0.05
% DG0 is the glitches selection for threshold = 0.1
if icorrect == 0
    %dirALL=[dirDROP 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGFP0_var1'];
    dirALL=[dirDROP 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGFP0_var1'];
    dirSEED=[dirDROP 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGSP0'];
    
else
    
    if istep == 1     %one step of deglitching
        dirALL=[dirDROP 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGF1'];
    elseif istep == 2 %two steps of deglitching
        dirALL=[dirDROP 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGF2'];      
    end
end

fc=0.005;
sfc=num2str(fc);
dt=.5;                     % output frequency of the signal
%%%%%%% attention, dans certains cas ( par exemple iday = 421), 
%%%%%%% il y a des decallage de plus de 10 msec a expliquer.
%%%%%%% on reduit ici la tolerance à 10% du samplin
%dt_tol=10e-3;
dt_tol=dt/10;

%% Dip and Azimuth + rotation matrix

thetaVBBu=-29.28;
thetaVBBv=-29.33;
thetaVBBw=-29.61;
phiVBBu=135.11;
phiVBBv=15.04;
phiVBBw=254.96;
Geom= [ sind(thetaVBBu) cosd(thetaVBBu)*cosd(phiVBBu) cosd(thetaVBBu)*sind(phiVBBu);
        sind(thetaVBBv) cosd(thetaVBBv)*cosd(phiVBBv) cosd(thetaVBBv)*sind(phiVBBv);
        sind(thetaVBBw) cosd(thetaVBBw)*cosd(phiVBBw) cosd(thetaVBBw)*sind(phiVBBw)];
recomp=inv(Geom);



%% Special case you read event catalog

fid = fopen('NewListEvent.txt');

S=textscan(fid,'%4d-%2d-%2dT%2d:%2d:%2d S%4d%s %s %s %s %s %4d-%2d-%2dT%2d:%2d:%2d %4d-%2d-%2dT%2d:%2d:%2d %s');
iyre=S{1};
imoe=S{2};
idae=S{3};
ihrz=S{4};
imne=S{5};
isee=S{6};
isoe=S{7};
ilee=S{8};
type=S{9};
iyrs=S{13};
imos=S{14};
idas=S{15};
ihrs=S{16};
imns=S{17};
ises=S{18};
iyrf=S{19};
imof=S{20};
idaf=S{21};
ihrf=S{22};
imnf=S{23};
isef=S{24};
Nevent=length(type);
%%

%% START DEGLITCHING PART

%%%%%% if iextract == 1 iev is the number of event
%%%%%%         else     iev is the sol number

% SOL 524 : veille du reglage TCDM (525)


for iev=solstart:solend                                                          % loop around the sols or the list event 
    if iextract == 1

        if iev ~= 601 & iev ~= 592 & iev ~= 609

            %%%%% Events types are
            %    'HIGHFR'
            %    'BROADB'
            %    'LOW_FR'
            %    '2.4_HZ'

            %%%%%%%%%% select the sol of the event
            if isequal(type(iev),{'BROADB'}) | isequal(type(iev),{'LOW_FR'})
                Tstart=datenum(double([iyrs(iev),imos(iev),idas(iev),ihrs(iev),imns(iev),ises(iev)]));
                Tend=datenum(double([iyrf(iev),imof(iev),idaf(iev),ihrf(iev),imnf(iev),isef(iev)]));
                iday=isoe(iev)
                sday=num2str(iday);
                Tstart=Tstart-1/24;
                Tend=Tend+1/24;
                extr=strcat('_S0',sday,sprintf('%s',ilee{iev}),'_');
                extr0=strcat('S0',sday,sprintf('%s',ilee{iev}));
                igo=1
            end
        end

    else
        igo = 1; % set parameters to selected sols
        iday=iev;
            sday=num2str(iday);
        extr='_sol_';
        extr0=strcat('sol',sday);
    end % end if iextract == 1

    if igo == 1

        if iday  >=74   & iday <= 167  % data were already @ 2sps
            dirU=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbu'];
            dirV=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbv'];
            dirW=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbw'];
        else iday >= 168               % data decimated from 20Hz to 2Hz 
            dirU=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbue'];
            dirV=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbve'];
            dirW=[dirDROP 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbwe'];
        end

    
    close all                                                                           % All this clears should be rearranged
    clearvars -except dirU dirV dirW dirALL dirDROP iday fc sfc recomp epsi dt dt_tol ...
        icorrect irot Geom Nhar ifilter ideriv istep valthreshold Tstart Tend extr iday iextract ...
        type isoe iyrs imos idas ihrs imns ises iyrf imof idaf ihrf imnf isef ilee iev extr0 dirSEED Nevent igo iseed iequal


    %% lecture des poles et zero 
    % Attention dans la calcul du gain, bien prendre uGain*uA0
    
    
    [ugain uzero upole units uA0]=read_resp_v2('RESP.XB.ELYSE.02.BHU_new');
    [vgain vzero vpole units vA0]=read_resp_v2('RESP.XB.ELYSE.02.BHV_new');
    [wgain wzero wpole units wA0]=read_resp_v2('RESP.XB.ELYSE.02.BHW_new');
    
   
    
    Uz=uzero;
    Up=upole;
    Vz=vzero;
    Vp=vpole;
    Wz=wzero;
    Wp=wpole;
    %%%%%%
    sday=num2str(iday);
    srot=num2str(irot);
    sthres=num2str(valthreshold);
    sstep=num2str(istep);
    sepsi=num2str(epsi);
    %% Read the data
    filedatau=[dirU,'/sol',sday,'.mat'];
    filedatav=[dirV,'/sol',sday,'.mat'];
    filedataw=[dirW,'/sol',sday,'.mat'];
    
    % Check sur l'existence
    if exist(filedatau) ~= 0 & exist(filedatav) ~= 0 & exist(filedataw) ~= 0
        % U
        load(filedatau);
        t1=t;ts1=ts;tlmst1=tlmst;
        % V
        load(filedatav);
        t2=t;ts2=ts;tlmst2=tlmst;
        % W
        load(filedataw);
        t3=t;ts3=ts;tlmst3=tlmst;

if length(ts1) ~= 0 & length(ts1) ~= 0 & length(ts1) ~= 0 

%%%%% check missing data and synchronize on common start of u,v,w

    kk1=1;
    while(ts1(kk1) == 0 | ts1(kk1) ~= ts1(kk1))
        kk1=kk1+1;
    end
    kk2=1;
    while(ts2(kk2) == 0 | ts2(kk2) ~= ts2(kk2)) 
        kk2=kk2+1;
    end
    kk3=1;
    while(ts3(kk3) == 0 | ts3(kk3) ~= ts3(kk3))
        kk3=kk3+1;
    end

    tstart=max([min(t1(kk1:end)),min(t2(kk2:end)),min(t3(kk3:end))]);

    nook=0;

    while(abs(t1(kk1)-tstart)*24*3600 > dt_tol)   %   Delta time smaller than 5 ms assume same time
        kk1=kk1+1;
        if kk1 > length(t1)
            nook=1;
            break;
        end
    end
    while(abs(t2(kk2)-tstart)*24*3600 > dt_tol)   %   Delta time smaller than 5 ms assume same time
        kk2=kk2+1;
        if kk2 > length(t2)
            nook=1;
            break;
        end
    end
    while(abs(t3(kk3)-tstart)*24*3600 > dt_tol)   %   Delta time smaller than 5 ms assume same time
        kk3=kk3+1;
        if kk3 > length(t3)
            nook=1;
            break;
        end
    end

if nook == 0

%%%%%%%%%%%%%%%%%%%
%%% detect any data interuption on u,v,w
Nu=length(t1); Nv=length(t2); Nw=length(t3);
for i=kk1+1:length(t1)
	if abs((t1(i)-t1(i-1))*24*3600 - dt) > dt_tol
		Nu=i-1;
		break;
	end
end
for i=kk2+1:length(t2)
	if abs((t2(i)-t2(i-1))*24*3600 - dt) > dt_tol
		Nv=i-1;
		break;
	end
end
for i=kk3+1:length(t3)
	if abs((t3(i)-t3(i-1))*24*3600 - dt) > dt_tol
		Nw=i-1;
		break;
	end
end
%%% take the shortest time
Tfin=[t1(Nu) t2(Nv) t3(Nw)];
[tz iz]=min(Tfin);
tfin=Tfin(iz);
%% define the end as the closest from shortest end
while( abs(t1(Nu)-tfin)*24*3600 > dt_tol) 
	Nu=Nu-1;
	if Nu < 1 
                nook=1;
                break;
        end
end
while( abs(t2(Nv)-tfin)*24*3600 > dt_tol) 
	Nv=Nv-1;
	if Nv < 1 
                nook=1;
                break;
        end
end
while( abs(t3(Nw)-tfin)*24*3600 > dt_tol) 
	Nw=Nw-1;
	if Nw < 1 
                nook=1;
                break;
        end
end

if nook == 0

%%%%%%%%%%%%%%%%%%%
utlmst=tlmst1(kk1:Nu);
ut=t1(kk1:Nu);
u=ts1(kk1:Nu);
vtlmst=tlmst2(kk2:Nv);
vt=t2(kk2:Nv);
v=ts2(kk2:Nv);
wtlmst=tlmst3(kk3:Nw);
wt=t3(kk3:Nw);
w=ts3(kk3:Nw);
%%%%%% raw data
%uraw=u;
%vraw=v;
%wraw=w;
%%%%% force interpolation on the longest components in case of holes
Nu=length(ut);
Nv=length(vt);
Nw=length(wt);
Ncase=[Nu, Nv, Nw];
[Nmax ic]=max(Ncase);
Ncase=Ncase-Nmax;
if max(abs(Ncase)) > 0
    warning1='ATTENTION LES COMPOSANTES U,V,W N ONT PAS LE MEME NOMBRE DE POINTS'
    if ic == 1
        v=interp1(vt,v,ut,'linear');
        vtlmst=interp1(vt,vtlmst,ut,'linear');
        vt=interp1(vt,vt,ut,'linear');
        w=interp1(wt,w,ut,'linear');
        wtlmst=interp1(wt,wtlmst,ut,'linear');
        wt=interp1(wt,wt,ut,'linear');
    elseif ic == 2
        u=interp1(ut,u,vt,'linear');
        utlmst=interp1(ut,utlmst,vt,'linear');
        ut=interp1(ut,ut,vt,'linear');
        w=interp1(wt,w,vt,'linear');
        wtlmst=interp1(wt,wtlmst,vt,'linear');
        wt=interp1(wt,wt,vt,'linear');
    elseif ic == 3
        v=interp1(vt,v,wt,'linear');
        vtlmst=interp1(vt,vtlmst,wt,'linear');
        vt=interp1(vt,vt,wt,'linear');
        u=interp1(ut,u,wt,'linear');
        utlmst=interp1(ut,utlmst,wt,'linear');
        ut=interp1(ut,ut,wt,'linear');
    end
end
%%%%%% flag start and end
urawstart=u(1);urawend=u(end);
vrawstart=v(1);vrawend=v(end);
wrawstart=w(1);wrawend=w(end);


%% Detrend du signal par une méthode de moindre carré sur la transformée de Fourier. 
% Nhar = 12 ce qui correspond (0.1mHz environ)
[TTu,bbTTu]=FourierleastsquarewithNaN(utlmst,u,Nhar);
for i=1:length(u)
    if u(i) == u(i)
        u(i)=u(i)-TTu(i);
    else
        u(i)=0.;
    end
end
%%%% v component
[TTv,bbTTv]=FourierleastsquarewithNaN(vtlmst,v,Nhar);
for i=1:length(v)
    if v(i) == v(i)
        v(i)=v(i)-TTv(i);
    else
        v(i)=0.;
    end
end
%%%% w component
[TTw,bbTTw]=FourierleastsquarewithNaN(wtlmst,w,Nhar);
for i=1:length(w)
    if w(i) == w(i)
        w(i)=w(i)-TTw(i);
    else
        w(i)=0.;
    end
end
%%%%% equalize the gain ( Here NaN must be removed)
if iequal == 1
    equalize
end
%%%%%%%%
if ideriv == 1
    %%%%%% jerk
    N=length(ut);
    NF=N/2+1;
    df=1./(N*dt);
    freq=[0:NF-1]'*df;
    clear spec;
    spec=fft(u);
    spec(1:NF)=spec(1:NF).*complex(0.,2*pi*freq(1:NF));
    du=ifft(spec,'symmetric');
    N=length(vt);
    NF=N/2+1;
    df=1./(N*dt);
    freq=[0:NF-1]'*df;
    clear spec;
    spec=fft(v);
    spec(1:NF)=spec(1:NF).*complex(0.,2*pi*freq(1:NF));
    dv=ifft(spec,'symmetric');
    N=length(wt);
    NF=N/2+1;
    df=1./(N*dt);
    freq=[0:NF-1]'*df;
    clear spec;
    spec=fft(w);
    spec(1:NF)=spec(1:NF).*complex(0.,2*pi*freq(1:NF));
    dw=ifft(spec,'symmetric');
    u=du;
    v=dv;
    w=dw;
end
%%%% filter
if ifilter == 2
    [B,A]=butter(4,.2*dt,'high');         %%%% high pass freq
    uf=filtfilt(B,A,u);            %%%%
    vf=filtfilt(B,A,v);            %%%%
    wf=filtfilt(B,A,w);            %%%%
    u=uf;
    v=vf;
    w=wf;
end
%%%%


if iextract == 1

    %%%%% take only the event

    [a iau]=min(abs(ut-Tstart));
    [b ibu]=min(abs(ut-Tend));
    [a iav]=min(abs(vt-Tstart));
    [b ibv]=min(abs(vt-Tend));
    [a iaw]=min(abs(wt-Tstart));
    [b ibw]=min(abs(wt-Tend));
    [ut,utlmst,u]=extractanderase0(iau,ibu,ut,utlmst,u);
    [vt,vtlmst,v]=extractanderase0(iav,ibv,vt,vtlmst,v);
    [wt,wtlmst,w]=extractanderase0(iaw,ibw,wt,wtlmst,w);

    %%%%%

else

    %%%%% take all the sol

    Tstart=min([ut(1),vt(1),wt(1)]);
    Tend=max([ut(end),vt(end),wt(end)]);

end

%%%% generate Z axis
uuz=recomp(1,1)*u+recomp(1,2)*v+recomp(1,3)*w;
uux=recomp(2,1)*u+recomp(2,2)*v+recomp(2,3)*w;
uuy=recomp(3,1)*u+recomp(3,2)*v+recomp(3,3)*w;
TTuuz=recomp(1,1)*TTu+recomp(1,2)*TTv+recomp(1,3)*TTw;
TTuux=recomp(2,1)*TTu+recomp(2,2)*TTv+recomp(2,3)*TTw;
TTuuy=recomp(3,1)*TTu+recomp(3,2)*TTv+recomp(3,3)*TTw;
bbTTuuz=recomp(1,1)*bbTTu+recomp(1,2)*bbTTv+recomp(1,3)*bbTTw;
bbTTuux=recomp(2,1)*bbTTu+recomp(2,2)*bbTTv+recomp(2,3)*bbTTw;
bbTTuuy=recomp(3,1)*bbTTu+recomp(3,2)*bbTTv+recomp(3,3)*bbTTw;
%%%%%%%%%%%%%%%%%%
if irot == 1
	u_old=u;
	v_old=v;
	w_old=w;
	u=uuz;
	v=uux;
	w=uuy;
end
%%%%%%%%%%%%%%%%%%

subplot(1,4,4)
plot(ut,uuz,'k')
%%%%
%NS=1024*5;
NS=1024;
NZERO=10000;
%%% .3 and .5 for full
%%% .04 for most
%%% .001 for calibration
dtthreshold=3*dt;
tfthreshold=valthreshold;
THRESHOLD3D=valthreshold;
igg(1:3)=0;
dgmin=2^25;
dgmax=-2^25;

udg=u;    % Duplication initiale des axes pour les envoyés dnas les fonctions
vdg=v;    % deglitch5 et deglitch7
wdg=w;
uuzdg=uuz;

% === sortie ===
%%%% si irot = 0, U,V,W sont U,V,W, si irot=1, U,V,W sont Z, N, E
%   unew:     : signal des glitchs en DU ( avec 0 si pas de glitch)
%   uglitch_1 : signal des glitchs en DU ( avec NaN si pas de glitch)    
%   udrift_1  : drift inverse mais pas retire
%   udg        :  signal degltiché
%   u0          : signal brute non degltiche mais sans tendance
%                        idem pour v, w
%   spredictor(:,1:3) 
%   uuznew, uuzdrift  pas vraiment utilisé
%   NGLITCH_1(:,1:3)  nombre de glitch 
%   Var0_1(:,1:3) variance residuelle pour tous les max
%   Var1_1 : rapport signal a glitch
%   val_1    : variance residulelle des glitch selectionne
%   val0_1  : variance residulelle de tous les glitchs detectés jusqu'a 50% de rediction
%   Amp_1 : amplitude du glitch Amp_1(:,numero de glitch, axe)
%   1ere dimention = NG2, nb de parametre du glitch ( reglé dans deglitch7)
%         deux cas: NG2= 5 ( 3 parametres de glitch + moyenne + tendance )
%   Delay_1: centre du glitch ( NG2=5) ou delai du glitch (NG2=3)
%   Green : reponse impulsionnelle utilisée
%   ND duree des fonctions de green
%   im position du max sur 1024 pt 
%   g0 
%   ND2: nombre de points avant le debut du glitch  ( longueur du glich via length(Green(:,1))
%   NDSUM et NDS variable cachée 
%   NZEROALL_A noombre de glitchs testes < 0.5
%   AmpSpike amplitude des spikes
%   jdx_1 indice des spike
% === entree ===
%   ut, utlmst, temps en UTC et LMST
%   u donnees
%   udg donnees entree
%   Up,Uz, pole et zero


[unew,uglitch_1,udrift_1,udiff_1,udg,u0,... 
    vnew,vglitch_1,vdrift_1, vdiff_1,vdg,v0,...
    wnew,wglitch_1,wdrift_1,wdiff_1,wdg,w0,...
    uuznew,uuzdrift_1,NGLITCH_1,Var0_1,Var1_1,val_1,val0_1,idx_1,...
    Amp_1,Delay_1,Green,ND,im,g0,ND2,NG1,NG2,...
    dgmin,dgmax,NZEROALL_1,spredictor_1,jdx_1,AmpSpike_1,reducA_1, Var]...
    =deglitch5_full_3axis(ut,utlmst,udg,u,...
    vt,vtlmst,vdg,v,...
    wt,wtlmst,wdg,w,...
    uuzdg,uuz,recomp,epsi,NS,dt,Uz,Up,Vz,Vp,Wz,Wp,NZERO,...
    valthreshold,dtthreshold,tfthreshold,dgmin,dgmax,1,iday,...
    icorrect,ideriv,igg,iequal);

%%%%% if no glitch detected, likely an error and jup


if max(NGLITCH_1) > 0

    uuzdg=recomp(1,1)*udg+recomp(1,2)*vdg+recomp(1,3)*wdg;

    if max(abs(udg)) > 0;
        dessine_spectro(utlmst,1000,1,1,dt,100,u,udg,'B')
    end
    if max(abs(vdg)) > 0;
        dessine_spectro(vtlmst,1000,1,2,dt,100,v,vdg,'B')
    end
    if max(abs(wdg)) > 0;
        dessine_spectro(wtlmst,1000,1,3,dt,100,w,wdg,'B')
    end

    udgdp=udg-spredictor_1(:,1);
    vdgdp=vdg-spredictor_1(:,2);
    wdgdp=wdg-spredictor_1(:,3);

    spredictor=spredictor_1;

    udg1=udg;
    vdg1=vdg;
    wdg1=wdg;

    if max(abs(udg)) > 0;
        dessine_spectro(utlmst,3000,1,1,dt,100,udg,udgdp,'B')
    end
    if max(abs(vdg)) > 0;
        dessine_spectro(vtlmst,3000,1,2,dt,100,vdg,vdgdp,'B')
    end
    if max(abs(wdg)) > 0;
        dessine_spectro(wtlmst,3000,1,3,dt,100,wdg,wdgdp,'B')
    end

    udg=udgdp;
    vdg=vdgdp;
    wdg=wdgdp;

    udg1=udg;
    vdg1=vdg;
    wdg1=wdg;

    if istep == 2                                                                    % 2ème passage de l'algo de détection (un glitch peut en cacher un autre)

        igg=NGLITCH_1;

        [unew,uglitch_2,udrift_2,udiff_2,udg,u0,...
            vnew,vglitch_2,vdrift_2,vdiff_2,vdg,v0,...
            wnew,wglitch_2,wdrift_2,wdiff_2,wdg,w0,...
            uuznew,uuzdrift_1,NGLITCH_2,Var0_2,Var1_2,val_2,val0_2,...
            idx_2,Amp_2,Delay_2,Green,ND,im,g0,ND2,NG1,NG2,dgmin,dgmax,...
            NZEROALL_2,spredictor_2,jdx_2,AmpSpike_2,reducA_2, Var]...
            =deglitch5_full_3axis(ut,utlmst,udg,udg,...
            vt,vtlmst,vdg,vdg,wt,wtlmst,wdg,wdg,...
            uuzdg,uuz,recomp,epsi,NS,dt,Uz,Up,Vz,Vp,Wz,Wp,...
            NZERO,valthreshold,dtthreshold,tfthreshold,...
            dgmin,dgmax,2,iday,icorrect,ideriv,igg,iequal);

    if abs(udg) > 0
    dessine_spectro(utlmst,1000,2,1,dt,100,udg1,udg,'B')
    end
    if abs(vdg) > 0
    dessine_spectro(vtlmst,1000,2,2,dt,100,vdg1,vdg,'B')
    end
    if abs(wdg) > 0
    dessine_spectro(wtlmst,1000,2,3,dt,100,wdg1,wdg,'B')
    end

    udgdp=udg-spredictor_2(:,1);
    vdgdp=vdg-spredictor_2(:,2);
    wdgdp=wdg-spredictor_2(:,3);

    spredictor=spredictor+spredictor_2;

    if abs(udg) > 0
        dessine_spectro(utlmst,3000,2,1,dt,300,udg,udgdp,'B')
    end
    if abs(vdg) > 0
        dessine_spectro(vtlmst,3000,2,2,dt,300,vdg,vdgdp,'B')
    end
    if abs(wdg) > 0
        dessine_spectro(wtlmst,3000,2,3,dt,300,wdg,wdgdp,'B')
    end

    NGLITCH=NGLITCH_1+NGLITCH_2;
    NZEROALL=NZEROALL_1+NZEROALL_2;

    udg2=udg;
    vdg2=vdg;
    wdg2=wdg;

    udg=udgdp;
    vdg=vdgdp;
    wdg=wdgdp;

    for iaxe=1:3

        if NGLITCH_1(iaxe) > 0
            val0(1:NGLITCH_1(iaxe),iaxe)=val0_1(1:NGLITCH_1(iaxe),iaxe);
            val(1:NGLITCH_1(iaxe),iaxe)=val_1(1:NGLITCH_1(iaxe),iaxe);
            idx(1:NGLITCH_1(iaxe),iaxe)=idx_1(1:NGLITCH_1(iaxe),iaxe);
            Amp(1:NG2,1:NGLITCH_1(iaxe),iaxe)=Amp_1(1:NG2,1:NGLITCH_1(iaxe),iaxe);
            Delay(1:NGLITCH_1(iaxe),iaxe)=Delay_1(1:NGLITCH_1(iaxe),iaxe);
            AmpSpike(1:NGLITCH_1(iaxe),iaxe)=AmpSpike_1(1:NGLITCH_1(iaxe),iaxe);
            reducA(1:NGLITCH_1(iaxe),iaxe)=reducA_1(1:NGLITCH_1(iaxe),iaxe);
        else
            val0(1,iaxe)=0.;
            val(1,iaxe)=0.;
            idx(1,iaxe)=0.;
            Amp(1,1,iaxe)=0.;
            ampSpike(1,iaxe)=0.;
            reducA(1,iaxe)=0.;
         end

        if NGLITCH_2(iaxe) > 0
            val0(NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=val0_2(1:NGLITCH_2(iaxe),iaxe);
            val(NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=val_2(1:NGLITCH_2(iaxe),iaxe);
            idx(NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=idx_2(1:NGLITCH_2(iaxe),iaxe);
            Amp(1:NG2,NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=Amp_2(1:NG2,1:NGLITCH_2(iaxe),iaxe);
            Delay(NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=Delay_2(1:NGLITCH_2(iaxe),iaxe);
            AmpSpike(NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=AmpSpike_2(1:NGLITCH_2(iaxe),iaxe);
            reducA(NGLITCH_1(iaxe)+1:NGLITCH_1(iaxe)+NGLITCH_2(iaxe),iaxe)=reducA_2(1:NGLITCH_2(iaxe),iaxe);
        end


    end

    if NGLITCH_1(1) > 0
        ujdx(1:NGLITCH_1(1))=jdx_1(1:NGLITCH_1(1),1)+ND2;
    else
        ujdx(1)=0;
    end
    if NGLITCH_1(2) > 0
        vjdx(1:NGLITCH_1(2))=jdx_1(1:NGLITCH_1(2),2)+ND2;
    else
        vjdx(1)=0;
    end
    if NGLITCH_1(3) > 0
        wjdx(1:NGLITCH_1(3))=jdx_1(1:NGLITCH_1(3),3)+ND2;
    else
        wjdx(1)=0;
    end
    
    if NGLITCH_2(1) > 0
        ujdx(NGLITCH_1(1)+1:NGLITCH_1(1)+NGLITCH_2(1))=jdx_2(1:NGLITCH_2(1),1)+ND2;
    end
    if NGLITCH_2(2) > 0
        vjdx(NGLITCH_1(2)+1:NGLITCH_1(2)+NGLITCH_2(2))=jdx_2(1:NGLITCH_2(2),2)+ND2;
    end
    if NGLITCH_2(3) > 0
        wjdx(NGLITCH_1(3)+1:NGLITCH_1(3)+NGLITCH_2(3))=jdx_2(1:NGLITCH_2(3),3)+ND2;
    end

    elseif istep == 1

        NGLITCH=NGLITCH_1;
        NZEROALL=NZEROALL_1;
        val0=val0_1;
        val=val_1;
        idx=idx_1;
        Amp=Amp_1;

        ujdx(1:NGLITCH_1(1))=jdx_1(1:NGLITCH_1(1),1)+ND2;
        vjdx(1:NGLITCH_1(2))=jdx_1(1:NGLITCH_1(2),2)+ND2;
        wjdx(1:NGLITCH_1(3))=jdx_1(1:NGLITCH_1(3),3)+ND2;
        AmpSpike=AmpSpike_1;
        reducA=reducA_1;

    end % end istep == 2 (2nd passage de l'algo)

    %%%%%% 3D analysis
    timeu=[0:ND-1]*dt;
    timev=[0:ND-1]*dt;
    timew=[0:ND-1]*dt;
    uNGLITCH=NGLITCH(1);
    vNGLITCH=NGLITCH(2);
    wNGLITCH=NGLITCH(3);
    uval(:)=val(:,1);
    uuval(:)=sort(val0(1:NGLITCH(1),1));
    vval(:)=val(:,2);
    vvval(:)=sort(val0(1:NGLITCH(2),2));
    wval(:)=val(:,3);
    wwval(:)=sort(val(1:NGLITCH(3),3));
    uAmp(:,:)=Amp(:,:,1);
    vAmp(:,:)=Amp(:,:,2);
    wAmp(:,:)=Amp(:,:,3);
    uidx(:)=idx(:,1);
    vidx(:)=idx(:,2);
    widx(:)=idx(:,3);
    uGreen(:,:)=Green(:,:,1);
    vGreen(:,:)=Green(:,:,2);
    wGreen(:,:)=Green(:,:,3);

    %%%% common part
    ztstart=max([ut(1),vt(1),wt(1)]);
    ztend=min([ut(end),vt(end),wt(end)]);
    [usdz isuz]=min(abs(ut-ztstart));
    [vsdz isvz]=min(abs(vt-ztstart));
    [wsdz iswz]=min(abs(wt-ztstart));
    [uedz ieuz]=min(abs(ut-ztend));
    [vedz ievz]=min(abs(vt-ztend));
    [wedz iewz]=min(abs(wt-ztend));
    NDz=max([ieuz-isuz+1 , ievz-isvz+1 , iewz-iswz+1]); 
    %%%%%% these are the outputs
    uzdg(1:NDz)=udg(isuz:ieuz);
    vzdg(1:NDz)=vdg(isvz:ievz);
    wzdg(1:NDz)=wdg(iswz:iewz);
    uzz(1:NDz)=u(isuz:ieuz);
    vzz(1:NDz)=v(isvz:ievz);
    wzz(1:NDz)=w(iswz:iewz);
    zt(1:NDz)=ut(isuz:ieuz);
    ztlmst(1:NDz)=utlmst(isuz:ieuz);
    %%%%%%%%%%%%%%%%%%%%%%%%ùù
    zflag(1:NDz)=0;
    varbest(1:NDz)=10.;
    UU(1:NDz)=0.;
    VV(1:NDz)=0.;
    WW(1:NDz)=0.;
    for i=1:uNGLITCH
    ii=uidx(i)-isuz+1;
    if ii > 0 & ii <= NDz
    varbest(ii)=min(varbest(ii),uval(i));
    UU(ii)=uAmp(1,i);
    end
    end
    for i=1:vNGLITCH
    ii=vidx(i)-isvz+1;
    if ii > 0 & ii <= NDz
    varbest(ii)=min(varbest(ii),vval(i));
    VV(ii)=vAmp(1,i);
    end
    end
    for i=1:wNGLITCH
    ii=widx(i)-iswz+1; 
    if ii > 0 & ii <= NDz
        varbest(ii)=min(varbest(ii),wval(i));
        WW(ii)=wAmp(1,i);
    end
    end

    valbest=sort(varbest);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    UUraw=UU;
    VVraw=VV;
    WWraw=WW;

    for i=1:NDz
        if varbest(i) < THRESHOLD3D
            zflag(i)=1;
            AMP=sqrt(UU(i)^2+VV(i)^2+WW(i)^2);
            if AMP == 0.
                AMP=1.;
            end
            UU(i)=UU(i)/AMP;
            VV(i)=VV(i)/AMP;
            WW(i)=WW(i)/AMP;
            %%%%% flag the weakly and non detected
            if abs(UU(i)) < 0.05 & UU(i) ~= 0.
                UU(i) = sign(UU(i)) * 0.05;
            end
            if abs(VV(i)) < 0.05 & VV(i) ~= 0.
                VV(i) = sign(VV(i)) * 0.05;
            end
            if abs(WW(i)) < 0.05 & WW(i) ~= 0.
                WW(i) = sign(WW(i)) * 0.05;
            end
        else
            varbest(i)=NaN(1);
            UU(i)=NaN(1);
            VV(i)=NaN(1);
            WW(i)=NaN(1);
        end
    end

    %%%% flag the accepted glitch
    for i=1:uNGLITCH
        uflag(i)=0;
        ii=uidx(i)-isuz+1 ;
        if ii > 0 & ii <= NDz & zflag(ii) == 1
            uflag(i)=1;
        else
            uval(i)=NaN(1);
        end
    end
    for i=1:vNGLITCH
        vflag(i)=0;
        ii=vidx(i)-isvz+1 ;
        if ii > 0 & ii <= NDz & zflag(ii) == 1
            vflag(i)=1;
        else
            vval(i)=NaN(1);
        end
    end
    
    for i=1:wNGLITCH
    wflag(i)=0;
    ii=widx(i)-iswz+1;
    if ii > 0 & ii <= NDz & zflag(ii) == 1
    wflag(i)=1;
    else
    wval(i)=NaN(1);
    end
    end

    %%%%%% forced removal with other axis informations

    ND=length(uGreen(:,1));

    unew_3(1:NDz)=0.;
    uglitch_3(1:NDz)=0.;
    udrift_3(1:NDz)=0.;
    for i=1:uNGLITCH
        if uflag(i) == 1
                for j=1:ND
                        for j1=1:NG1
                                unew_3(uidx(i)+j-isuz)=unew_3(uidx(i)+j-isuz)+uAmp(j1,i)*uGreen(j,j1);
                        end
                        for j1=1:NG2
                                uglitch_3(uidx(i)+j-isuz)=uglitch_3(uidx(i)+j-isuz)+uAmp(j1,i)*uGreen(j,j1);
                        end
                        for j1=NG1+1:NG2
                                udrift_3(uidx(i)+j-isuz)=udrift_3(uidx(i)+j-isuz)+uAmp(j1,i)*uGreen(j,j1);
                        end
                end
        end
    end

    un=u(isuz:ieuz);

    %%% integrate the glitch
    clear spec
    spec=fft(unew_3);
    NF=floor(length(unew_3)/2)+1;
    df=1./dt/length(unew_3);
    ff=[0:NF-1]*df;
    for i=1:NF
        spec(i)=1./complex(0.,ff(i)*2*pi)*spec(i);
    end
    %spec(1)=0.;
    unew_3=ifft(spec,'symmetric');
    unn=un-unew_3';

    figure(500)
    plot(un,'k')
    hold on
    plot(unn,'r')

    unew_3=unew_3';
    uglitch_3=uglitch_3';

    %%%%%%%%%%

    figure(10)
    subplot(1,3,1)
    plot(datetime(datevec(ut)),u,'k')
    hold on
    plot(datetime(datevec(zt)),unn,'g')
    plot(datetime(datevec(ut)),udg,'k')
    plot([datetime(datevec(ut(1))),datetime(datevec(ut(end)))],[udg(1),udg(1)],'r','LineWidth',1)
    subplot(1,3,2)
    plot(datetime(datevec(vt)),v,'b')
    hold on
    plot(datetime(datevec(vt)),vdg,'b')
    plot([datetime(datevec(vt(1))),datetime(datevec(vt(end)))],[vdg(1),vdg(1)],'r','LineWidth',1)
    subplot(1,3,3)
    plot(datetime(datevec(wt)),w,'g')
    hold on
    plot(datetime(datevec(wt)),wdg,'g')
    plot([datetime(datevec(wt(1))),datetime(datevec(wt(end)))],[wdg(1),wdg(1)],'r','LineWidth',1)

 

    figure(32)
    if uNGLITCH > 0
        subplot(3,3,1)
        plot(datetime(datevec(ut(uidx(1:uNGLITCH)))),uAmp(1,1:uNGLITCH),'k*')
        subplot(3,3,2)
        plot(datetime(datevec(ut(uidx(1:uNGLITCH)))),uAmp(2,1:uNGLITCH)./uAmp(1,1:uNGLITCH),'k*')
        hold on
        plot([datetime(datevec(ut(1))),datetime(datevec(ut(end)))],[-.5, -.5],'-r')
        plot([datetime(datevec(ut(1))),datetime(datevec(ut(end)))],[.5, .5],'-r')
        %xlim([datetime(datevec(ut(1))),datetime(datevec(ut(end)))])
        ylim([-1 1])
    end
    if vNGLITCH > 0
        subplot(3,3,3)
        plot(datetime(datevec(vt(vidx(1:vNGLITCH)))),vAmp(3,1:vNGLITCH)./vAmp(1,1:vNGLITCH),'k*')
        %%%%%
        subplot(3,3,4)
        plot(datetime(datevec(vt(vidx(1:vNGLITCH)))),vAmp(1,1:vNGLITCH),'k*')
        subplot(3,3,5)
        plot(datetime(datevec(vt(vidx(1:vNGLITCH)))),vAmp(2,1:vNGLITCH)./vAmp(1,1:vNGLITCH),'k*')
        hold on
        plot([datetime(datevec(vt(1))),datetime(datevec(vt(end)))],[-.5, -.5],'-r')
        plot([datetime(datevec(vt(1))),datetime(datevec(vt(end)))],[.5, .5],'-r')
        %xlim([datetime(datevec(vt(1))),datetime(datevec(vt(end)))])
        ylim([-1 1])
    end
    if wNGLITCH > 0
        subplot(3,3,6)
        plot(datetime(datevec(wt(widx(1:wNGLITCH)))),wAmp(3,1:wNGLITCH)./wAmp(1,1:wNGLITCH),'k*')
        subplot(3,3,7)
        plot(datetime(datevec(wt(widx(1:wNGLITCH)))),wAmp(1,1:wNGLITCH),'k*')
        subplot(3,3,8)
        plot(datetime(datevec(wt(widx(1:wNGLITCH)))),wAmp(2,1:wNGLITCH)./wAmp(1,1:wNGLITCH),'k*')
        hold on
        plot([datetime(datevec(wt(1))),datetime(datevec(wt(end)))],[-.5, -.5],'-r')
        plot([datetime(datevec(wt(1))),datetime(datevec(wt(end)))],[.5, .5],'-r')
        %xlim([datetime(datevec(wt(1))),datetime(datevec(wt(end)))])
        ylim([-1 1])
        subplot(3,3,9)
        plot(datetime(datevec(wt(widx(1:wNGLITCH)))),wAmp(3,1:wNGLITCH)./wAmp(1,1:wNGLITCH),'k*')
        %%%%%%
        xlim([datetime(datevec(wt(1))),datetime(datevec(wt(end)))])
    end

    %%%%%%% 

    % Z(nadir), X - North, Y -East
    if irot == 1
        ZZ=UU;
        XX=VV;
        YY=WW;
        UU=Geom(1,1)*ZZ+Geom(1,2)*XX+Geom(1,3)*YY;
        VV=Geom(2,1)*ZZ+Geom(2,2)*XX+Geom(2,3)*YY;
        WW=Geom(3,1)*ZZ+Geom(3,2)*XX+Geom(3,3)*YY;
    %%%%% rescale
        NORM=sqrt(UU.*UU+VV.*VV+WW.*WW);
        UU=UU./NORM;
        VV=VV./NORM;
        WW=WW./NORM;
        ZZ=ZZ/(ugain*uA0);
        XX=XX/(ugain*uA0);
        YY=YY/(ugain*uA0);
    else
        ZZ=recomp(1,1)*UU/(ugain*uA0)+recomp(1,2)*VV/(ugain*uA0)+recomp(1,3)*WW/(ugain*uA0);
        XX=recomp(2,1)*UU/(ugain*uA0)+recomp(2,2)*VV/(ugain*uA0)+recomp(2,3)*WW/(ugain*uA0);
        YY=recomp(3,1)*UU/(ugain*uA0)+recomp(3,2)*VV/(ugain*uA0)+recomp(3,3)*WW/(ugain*uA0);
    end

    iproject = 1;

    if iproject == 1

    iu=0;iv=0;iw=0;io=0;
    UUui(1)=0.;
    UUoi(1)=0.;
    UUc(1)=0.;
    for i=2:length(UUraw)
            UUc(i)=UUc(i-1)+uzdg(i);
            if abs(UU(i)) > .9 
                UUui(i) = UUui(i-1)+UUraw(i);
                UUoi(i) = UUoi(i-1);
            else
                UUui(i) = UUui(i-1);
                UUoi(i) = UUoi(i-1)+UUraw(i);
            end
    end
    %%%%%%%
    ztu=[];ztlmstu=[];ZZu=[];YYu=[];XXu=[];DIPu=[];AZIu=[];
    ztv=[];ztlmstv=[];ZZv=[];YYv=[];XXv=[];DIPv=[];AZIv=[];
    ztw=[];ztlmstw=[];ZZw=[];YYw=[];XXw=[];DIPw=[];AZIw=[];
    zto=[];ztlmsto=[];ZZo=[];YYo=[];XXo=[];DIPo=[];AZIo=[];
    for i=1:length(XX)
        norm=ZZ(i)*ZZ(i)+XX(i)*XX(i)+YY(i)*YY(i);
        if norm ~= 0
    %%%%%
            if abs(UU(i)) > .9
            norm=sqrt(norm);
            iu=iu+1;
            ztu(iu)=zt(i);
            ztlmstu(iu)=ztlmst(i);
            ZZu(iu)=ZZ(i)/norm;
            XXu(iu)=XX(i)/norm;
            YYu(iu)=YY(i)/norm;
            DIPu(iu)=atan2d(ZZ(i),sqrt(XX(i)*XX(i)+YY(i)*YY(i)));
            azii=atan2d(XX(i),YY(i));
            if azii < 0
                AZIu(iu,1) = NaN(1);
                AZIu(iu,2) = azii+180;
            else
                AZIu(iu,1) = azii;
                AZIu(iu,2) = NaN(1);
            end
            elseif abs(VV(i)) > .9
            norm=sqrt(norm);
            iv=iv+1;
            ztv(iv)=zt(i);
            ztlmstv(iv)=ztlmst(i);
            ZZv(iv)=ZZ(i)/norm;
            XXv(iv)=XX(i)/norm;
            YYv(iv)=YY(i)/norm;
            DIPv(iv)=atan2d(ZZ(i),sqrt(XX(i)*XX(i)+YY(i)*YY(i)));
            azii=atan2d(XX(i),YY(i));
            if azii < 0
                AZIv(iv,2) = azii+180;
                AZIv(iv,1) = NaN(1);
            else
                AZIv(iv,1) = azii;
                AZIv(iv,2) = NaN(1);
            end
            elseif abs(WW(i)) > .9
            norm=sqrt(norm);
            iw=iw+1;
            ztw(iw)=zt(i);
            ztlmstw(iw)=ztlmst(i);
            ZZw(iw)=ZZ(i)/norm;
            XXw(iw)=XX(i)/norm;
            YYw(iw)=YY(i)/norm;
            DIPw(iw)=atan2d(ZZ(i),sqrt(XX(i)*XX(i)+YY(i)*YY(i)));
            azii=atan2d(XX(i),YY(i));
            if azii < 0
                AZIw(iw,2) = azii+180;
                AZIw(iw,1) = NaN(1);
            else
                AZIw(iw,1) = azii;
                AZIw(iw,2) = NaN(1);
            end
            elseif norm == norm
            norm=sqrt(norm);
            io=io+1;
            zto(io)=zt(i);
            ztlmsto(io)=ztlmst(i);
            ZZo(io)=ZZ(i)/norm;
            XXo(io)=XX(i)/norm;
            YYo(io)=YY(i)/norm;
            DIPo(io)=atan2d(ZZ(i),sqrt(XX(i)*XX(i)+YY(i)*YY(i)));
            azii=atan2d(XX(i),YY(i));
            if azii < 0
                AZIo(io,2) = azii+180;
                AZIo(io,1) = NaN(1);
            else
                AZIo(io,1) = azii;
                AZIo(io,2) = NaN(1);
            end
            end
        end
    end
    figure(100*(istep-1)+35)
    if iu > 0
        subplot(1,3,1)
        plot(datetime(datevec(ztu)),DIPu,'r*')
        hold on
        subplot(1,3,2)
        plot(datetime(datevec(ztu)),AZIu(:,1),'r>')
        hold on
        plot(datetime(datevec(ztu)),AZIu(:,2),'r<')
        subplot(1,3,3)
        plot3(XXu,YYu,ZZu,'r*')
        hold on
    end
    if iv > 0
        subplot(1,3,1)
        plot(datetime(datevec(ztv)),DIPv,'g*')
        hold on
        subplot(1,3,2)
        plot(datetime(datevec(ztv)),AZIv(:,1),'g>')
        hold on
        plot(datetime(datevec(ztv)),AZIv(:,2),'g<')
        subplot(1,3,3)
        plot3(XXv,YYv,ZZv,'g*')
        hold on
    end
    if iw > 0
        subplot(1,3,1)
        plot(datetime(datevec(ztw)),DIPw,'b*')
        hold on
        subplot(1,3,2)
        plot(datetime(datevec(ztw)),AZIw(:,1),'b>')
        hold on
        plot(datetime(datevec(ztw)),AZIw(:,2),'b<')
        subplot(1,3,3)
        plot3(XXw,YYw,ZZw,'b*')
        hold on
    end
    subplot(1,3,1)
    plot(datetime(datevec(zto)),DIPo,'k*')
    hold on
    title('Dip angle')
    legend('U','V','W','Other')
    subplot(1,3,2)
    if length(zto) > 0
        plot(datetime(datevec(zto)),AZIo(:,1),'k>')
        hold on
        plot(datetime(datevec(zto)),AZIo(:,2),'k<')
        title('Azimuth angle')
        subplot(1,3,3)
        plot3(XXo,YYo,ZZo,'k*')
        hold on
    end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(100*(istep-1)+33)
    subplot(2,1,1)
    plot(datetime(datevec(zt)),UU,'k*')
    hold on
    plot(datetime(datevec(zt)),VV,'b*')
    plot(datetime(datevec(zt)),WW,'g*')
    %xlim([datetime(datevec(zt(1))),datetime(datevec(zt(end)))])
    subplot(2,1,2)
    semilogy(datetime(datevec(ut(uidx(1:uNGLITCH)))),uuval(1:uNGLITCH),'k*')
    hold on
    semilogy(datetime(datevec(vt(vidx(1:vNGLITCH)))),vvval(1:vNGLITCH),'b*')
    semilogy(datetime(datevec(wt(widx(1:wNGLITCH)))),wwval(1:wNGLITCH),'g*')
    semilogy(datetime(datevec(zt)),varbest,'mO')
    uVar0_limit(1:length(uval))=valthreshold;
    semilogy(datetime(datevec(ut(uidx(1:uNGLITCH)))),uVar0_limit(1:uNGLITCH),'r')
    uVar0_limit(1:length(uval))=THRESHOLD3D;
    semilogy(datetime(datevec(ut(uidx(1:uNGLITCH)))),uVar0_limit(1:uNGLITCH),'r')
    %xlim([datetime(datevec(zt(1))),datetime(datevec(zt(end)))])
    legend('Glitch U','Glitch V','Glitch W','Acceptation threshold','Accepted glitch')
    title('Glitch variance reduction')
    xlabel('Number of glitch tested')
    figure(100*(istep-1)+34)
    loglog(uuval,'r')
    hold on
    loglog(vvval,'g')
    loglog(wwval,'b')
    loglog(valbest,'k')
    ylim([1e-5 valthreshold]) 
    title('Variance reduction ranking')
    xlabel('Glitch number')
    ylabel('Variance after deglitching')
    legend('U','V','W','Best')

    %%%% final data after deglitching

    %%%% go back to right coordinates
    % udg,vdg,wdg  deglitched u,v,w in DU
    % u,v,w  raw u,v,w in DU
    % zdg,xdg,ydg  deglitched z,x,y in DU
    % uzdg,vzdg,wzdg  deglitched u,v,w in DU on the synchronized part with U TF
    % uzz,vzz,wzz  raw u,v,w in DU on the synchronized part with U TF
    % zzdg,xzdg,yzdg  deglitched z,x,y in DU on the synchronized part with U TF
    % zzz,xzz,yzz  raw z,x,y in DU on the synchronized part with U TF

    if irot == 1
        %%%%% correct the data as made in z,x,u
        zdg=udg;
        xdg=vdg;
        ydg=wdg;
        zpredictor=spredictor;
        zzz=uzz;
        xzz=vzz;
        yzz=wzz;
        zzdg=uzdg;
        xzdg=vzdg;
        yzdg=wzdg;

        if iseed == 1

            nameadd=['&',srot,'_',sthres,'_',sstep,'_',sepsi,'_',extr0,'_IPGP'];

            filename3=[dirSEED,'/XB.ELYSE.20.BHZ.',nameadd];
            filename4=[dirSEED,'/XB.ELYSE.21.BHZ.',nameadd];
            mkmseed_IPGP(filename3,zzdg,zt,2);
            mkmseed_IPGP(filename4,zzz,zt,2);

            filename3=[dirSEED,'/XB.ELYSE.20.BHN.',nameadd];
            filename4=[dirSEED,'/XB.ELYSE.21.BHN.',nameadd];
            mkmseed_IPGP(filename3,xzdg,zt,2);
            mkmseed_IPGP(filename4,xzz,zt,2);

            filename3=[dirSEED,'/XB.ELYSE.20.BHE.',nameadd];
            filename4=[dirSEED,'/XB.ELYSE.21.BHE.',nameadd];
            mkmseed_IPGP(filename3,yzdg,zt,2);
            mkmseed_IPGP(filename4,yzz,zt,2);

        end

        % generate u,v,w data fro z,x,y
        uzz=Geom(1,1)*zzz+Geom(1,2)*xzz+Geom(1,3)*yzz;
        vzz=Geom(2,1)*zzz+Geom(2,2)*xzz+Geom(2,3)*yzz;
        wzz=Geom(3,1)*zzz+Geom(3,2)*xzz+Geom(3,3)*yzz;
        uzdg=Geom(1,1)*zzdg+Geom(1,2)*xzdg+Geom(1,3)*yzdg;
        vzdg=Geom(2,1)*zzdg+Geom(2,2)*xzdg+Geom(2,3)*yzdg;
        wzdg=Geom(3,1)*zzdg+Geom(3,2)*xzdg+Geom(3,3)*yzdg;
        spredictor(:,1)=Geom(1,1)*zpredictor(:,1)+Geom(1,2)*zpredictor(:,2)+Geom(1,3)*zpredictor(:,3);
        spredictor(:,2)=Geom(2,1)*zpredictor(:,1)+Geom(2,2)*zpredictor(:,2)+Geom(2,3)*zpredictor(:,3);
        spredictor(:,3)=Geom(3,1)*zpredictor(:,1)+Geom(3,2)*zpredictor(:,2)+Geom(3,3)*zpredictor(:,3);
        udg=uzdg;
        vdg=vzdg;
        wdg=wzdg;
        %%%%
        u=u_old;
        v=v_old;
        w=w_old;
    else
        %%%% go back to right coordinates
        % udg,vdg,wdg  deglitched u,v,w in DU
        % u,v,w  raw u,v,w in DU
        % zdg,xdg,ydg  deglitched z,x,y in DU
        % uzdg,vzdg,wzdg  deglitched u,v,w in DU on the synchronized part with U TF
        % uzz,vzz,wzz  raw u,v,w in DU on the synchronized part with U TF
        % zzdg,xzdg,yzdg  deglitched z,x,y in DU on the synchronized part with U TF
        % zzz,xzz,yzz  raw z,x,y in DU on the synchronized part with U TF
        zdg=recomp(1,1)*udg+recomp(1,2)*vdg+recomp(1,3)*wdg;
        xdg=recomp(2,1)*udg+recomp(2,2)*vdg+recomp(2,3)*wdg;
        ydg=recomp(3,1)*udg+recomp(3,2)*vdg+recomp(3,3)*wdg;
        zpredictor(:,1)=recomp(1,1)*spredictor(:,1)+recomp(1,2)*spredictor(:,2)+recomp(1,3)*spredictor(:,3);
        zpredictor(:,2)=recomp(2,1)*spredictor(:,1)+recomp(2,2)*spredictor(:,2)+recomp(2,3)*spredictor(:,3);
        zpredictor(:,3)=recomp(3,1)*spredictor(:,1)+recomp(3,2)*spredictor(:,2)+recomp(3,3)*spredictor(:,3);
        zzdg=recomp(1,1)*uzdg+recomp(1,2)*vzdg+recomp(1,3)*wzdg;
        xzdg=recomp(2,1)*uzdg+recomp(2,2)*vzdg+recomp(2,3)*wzdg;
        yzdg=recomp(3,1)*uzdg+recomp(3,2)*vzdg+recomp(3,3)*wzdg;
        zzz=recomp(1,1)*uzz+recomp(1,2)*vzz+recomp(1,3)*wzz;
        xzz=recomp(2,1)*uzz+recomp(2,2)*vzz+recomp(2,3)*wzz;
        yzz=recomp(3,1)*uzz+recomp(3,2)*vzz+recomp(3,3)*wzz;
    end


    %udg=remove_equalizing(udg,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    vdg=remove_equalizing(vdg,Vz,Vp,vA0,vgain,Uz,Up,uA0,ugain,dt);
    wdg=remove_equalizing(wdg,Wz,Wp,wA0,wgain,Uz,Up,uA0,ugain,dt);
    %temp=remove_equalizing(spredictor(:,1),Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt); %spredictor(:,1)=temp;
    temp=remove_equalizing(spredictor(:,2),Vz,Vp,vA0,vgain,Uz,Up,uA0,ugain,dt); spredictor(:,2)=temp;
    temp=remove_equalizing(spredictor(:,3),Wz,Wp,wA0,wgain,Uz,Up,uA0,ugain,dt); spredictor(:,3)=temp;
    %uzdg=remove_equalizing(uzdg,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    vzdg=remove_equalizing(vzdg,Vz,Vp,vA0,vgain,Uz,Up,uA0,ugain,dt);
    wzdg=remove_equalizing(wzdg,Wz,Wp,wA0,wgain,Uz,Up,uA0,ugain,dt);
    %uzz=remove_equalizing(uzz,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    vzz=remove_equalizing(vzz,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    wzz=remove_equalizing(wzz,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    %zzdg=remove_equalizing(zzdg,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    %xzdg=remove_equalizing(xzdg,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    %yzdg=remove_equalizing(yzdg,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    %zzz=remove_equalizing(zzz,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    %xzz=remove_equalizing(xzz,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);
    %yzz=remove_equalizing(yzz,Uz,Up,uA0,ugain,Uz,Up,uA0,ugain,dt);

    %%% go in physical unit

    %%%% first remove the thermal signal 

    [uzzc] = remove_tfvitdt2(uzz',zt,dt,ugain,uzero,upole,uA0,fc);
    [vzzc] = remove_tfvitdt2(vzz',zt,dt,ugain,uzero,upole,uA0,fc);
    [wzzc] = remove_tfvitdt2(wzz',zt,dt,ugain,uzero,upole,uA0,fc);
    [uzdgc] = remove_tfvitdt2(uzdg',zt,dt,ugain,uzero,upole,uA0,fc);
    [vzdgc] = remove_tfvitdt2(vzdg',zt,dt,ugain,uzero,upole,uA0,fc);
    [wzdgc] = remove_tfvitdt2(wzdg',zt,dt,ugain,uzero,upole,uA0,fc);

    zz=recomp(1,1)*uzzc+recomp(1,2)*vzzc+recomp(1,3)*wzzc;
    xx=recomp(2,1)*uzzc+recomp(2,2)*vzzc+recomp(2,3)*wzzc;
    yy=recomp(3,1)*uzzc+recomp(3,2)*vzzc+recomp(3,3)*wzzc;
    zzdg=recomp(1,1)*uzdgc+recomp(1,2)*vzdgc+recomp(1,3)*wzdgc;
    xxdg=recomp(2,1)*uzdgc+recomp(2,2)*vzdgc+recomp(2,3)*wzdgc;
    yydg=recomp(3,1)*uzdgc+recomp(3,2)*vzdgc+recomp(3,3)*wzdgc;

    figure(200)
    subplot(3,1,1)
    plot(datetime(datevec(zt)),uzzc,'k')
    hold on
    plot(datetime(datevec(zt)),vzzc,'b')
    plot(datetime(datevec(zt)),wzzc,'g')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    title([sfc,'Hz HP data with glitch in acceleration'])
    ylabel('Acceleration m/s^2')
    subplot(3,1,2)
    plot(datetime(datevec(zt)),uzdgc,'k')
    hold on
    plot(datetime(datevec(zt)),vzdgc,'b')
    plot(datetime(datevec(zt)),wzdgc,'g')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    legend('U','V','W')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    title('UVW Deglitched in acceleration')
    ylabel('Acceleration m/s^2')
    subplot(3,1,3)
    %%%%%
    plot(datetime(datevec(zt)),zzdg,'r')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    hold on
    plot(datetime(datevec(zt)),-1.5e-8+xxdg,'b')
    plot(datetime(datevec(zt)),-3.e-8+yydg,'g')
    %legend('Z','X','Y','Glitch X','Glitch Y','Glitch Z');
    ylabel('Acceleration m/s^2')
    title('XYZ Deglitched in acceleration')

    figure(1300)

    subplot(4,1,1)
    plot(datetime(datevec(zt)),zzdg,'k')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    hold on
    plot(datetime(datevec(zt)),-1.5e-8+xxdg,'b')
    plot(datetime(datevec(zt)),-3.e-8+yydg,'g')
    legend('Z','X','Y');
    title('Deglitched - non filtered data')

    [B,A]=butter(5,.1*dt,'low');         %%%% high pass freq
    zs_c3=filtfilt(B,A,zzdg);            %%%%
    xs_c3=filtfilt(B,A,xxdg);            %%%%
    ys_c3=filtfilt(B,A,yydg);            %%%%

    [B,A]=butter(5,1/15*dt,'high');         %%%% high pass freq
    zs_c4=filtfilt(B,A,zs_c3);            %%%%
    xs_c4=filtfilt(B,A,xs_c3);            %%%%
    ys_c4=filtfilt(B,A,ys_c3);            %%%%

    subplot(4,1,2)
    plot(datetime(datevec(zt)),zs_c4,'k')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    hold on
    plot(datetime(datevec(zt)),-1.5e-9+xs_c4,'b')
    plot(datetime(datevec(zt)),-3.e-9+ys_c4,'g')
    legend('Z','X','Y');
    title('Deglitched - [10sec-20sec] filtered data')

    subplot(4,1,3)
    Tshift=ND2*dt/(24*3600);

    plot(datetime(datevec(zt+Tshift)),UU,'*k')
    hold on
    plot(datetime(datevec(zt+Tshift)),VV,'*b')
    plot(datetime(datevec(zt+Tshift)),WW,'*g')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    title('Glitch director cosine')

    Nu=length(UU);
    zmark(1:Nu)=NaN(1);
    xmark(1:Nu)=NaN(1);
    ymark(1:Nu)=NaN(1);
    zmark2(1:Nu)=NaN(1);
    xmark2(1:Nu)=NaN(1);
    ymark2(1:Nu)=NaN(1);

    for i=1:Nu
    if UU(i) == UU(i)
    zmark(i)=zs_c4(i);
    zmark2(i)=zzdg(i);
    end
    if VV(i) == VV(i)
    xmark(i)=xs_c4(i);
    xmark2(i)=xxdg(i);
    end
    if WW(i) == WW(i)
    umark(i)=ys_c4(i);
    ymark2(i)=yydg(i);
    end
    end

    %%%% glitch marks
    figure(200)
    subplot(3,1,3)
    plot(datetime(datevec(zt+Tshift)),zmark2,'*c')
    hold on
    plot(datetime(datevec(zt+Tshift)),-1.5e-8+xmark2,'*c')
    plot(datetime(datevec(zt+Tshift)),-3.e-8+ymark2,'*c')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])

    figure(1300)
    subplot(4,1,2)
    plot(datetime(datevec(zt+Tshift)),zmark,'*c')
    %xlim([datetime(datevec(Tstart)) datetime(datevec(Tend))])
    hold on
    plot(datetime(datevec(zt+Tshift)),-1.5e-9+xmark,'*c')
    plot(datetime(datevec(zt+Tshift)),-3.e-9+ymark,'*c')

    subplot(4,1,1)
    plot(datetime(datevec(zt+Tshift)),zmark2,'*c')
    hold on
    plot(datetime(datevec(zt+Tshift)),-1.5e-8+xmark2,'*c')
    plot(datetime(datevec(zt+Tshift)),-3.e-8+ymark2,'*c')

    %%%%%
    %%vertical
    dessine_spectro(ztlmst,1100,1,1,dt,300,zz,zzdg,'A')
    dessine_spectro(ztlmst,1100,1,2,dt,300,xx,xxdg,'A')
    dessine_spectro(ztlmst,1100,1,3,dt,300,yy,yydg,'A')

    %%%%%% output

    for iaxe=1:3
        for i=1:length(ztlmst)
            if spredictor(i,iaxe) == 0.
                spredictor(i,iaxe) = NaN;
                %spredictor(i,iaxe) = 0;
            end
        end
    end
    
    
    Nuu=length(u);
    for i=1:Nuu
        glitchu(i)=u(i)-udg(i)-spredictor(i,1);
    end
    
    
    Nvu=length(v);
        for i=1:Nvu
            glitchv(i)=v(i)-vdg(i)-spredictor(i,2);
        end
        
     Nwu=length(w);
        for i=1:Nwu
            glitchw(i)=w(i)-wdg(i)-spredictor(i,3);
        end  
    
    filesave=[dirALL,'/XYZ',srot,'_',sthres,'_',sstep,'_',sepsi,extr,sday,'.mat'];
    save(filesave,'zt','ztlmst','uzz','vzz','wzz','uzdg','vzdg','wzdg','zz','xx','yy','zzdg','xxdg','yydg',...
        'fc','Nhar','isuz','ieuz','isvz','ievz','isvz','ievz','TTuuz','TTuux','TTuuy',...
        'bbTTuuz','bbTTuux','bbTTuuy','spredictor','Geom','recomp')
    
    
    filesave=[dirALL,'/UVW',srot','_',sthres,'_',sstep,'_',sepsi,extr,sday,'.mat'];
    save(filesave,'u','v','w','uzz','vzz','wzz','udg','ut','vdg','vt','wdg','wt','uzzc','vzzc','wzzc',...
        'uzdgc','vzdgc','wzdgc','fc','utlmst','vtlmst','wtlmst','spredictor', 'unew', 'glitchu', 'glitchv', 'glitchw', 'Var')
    
    
    filesave=[dirALL,'/G',srot','_',sthres,'_',sstep,'_',sepsi,extr,sday,'.mat'];
    save(filesave,'UU','VV','WW','XX','YY','ZZ','Tshift','xmark','ymark','zmark',...
        'xmark2','ymark2','zmark2','ztu','ztv','ztw','zto',...
        'DIPu','DIPv','DIPw','DIPo','AZIu','AZIv','AZIw','AZIo',...
        'ZZu','ZZv','ZZw','XXu','XXv','XXw','YYu','YYv','YYw',...
        'uzzc','vzzc','wzzc','uzdgc','vzdgc','wzdgc','fc',...
    'ztlmstu','ztlmstv','ztlmstw','ztlmsto','XXo','YYo','XXo','ZZo',...
        'uval','vval','wval','valbest','varbest','uidx','vidx','widx',...
        'uNGLITCH','vNGLITCH','wNGLITCH','ujdx','vjdx','wjdx','AmpSpike',...
        'reducA','iu','iv','iw','io','Amp','Delay','im',...
        'ND2','NG1','NG2','ND','Green' )

    %%%%%% save in mseed U in DU
    if iseed == 1
        nameadd=['&',srot,'_',sthres,'_',sstep,'_',sepsi,'_',extr0,'_IPGP'];
        Nuu=length(u);
        for i=1:Nuu
            glitchu(i)=u(i)-udg(i)-spredictor(i,1);
        end
        uuval(1:Nuu)=1.;
        uuval(uidx(1:uNGLITCH))=uval(1:uNGLITCH);
        filename1=[dirSEED,'/XB.ELYSE.20.BHU',nameadd];
        filename2=[dirSEED,'/XB.ELYSE.21.BHU',nameadd];
        filename3=[dirSEED,'/XB.ELYSE.22.BHU',nameadd];
        filename4=[dirSEED,'/XB.ELYSE.23.BHU',nameadd];
        filename5=[dirSEED,'/XB.ELYSE.24.BHU',nameadd];
        mkmseed_IPGP(filename1,glitchu,ut,2,'onefile');
        mkmseed_IPGP(filename2,spredictor(:,1),ut,2,'onefile');
        mkmseed_IPGP(filename3,udg',ut,2,'onefile');
        mkmseed_IPGP(filename4,u',ut,2,'onefile');
        mkmseed_IPGP(filename5,uuval',ut,2,'onefile');

        %%%%%% save in mseed V in DU
        Nvu=length(v);
        for i=1:Nvu
            glitchv(i)=v(i)-vdg(i)-spredictor(i,2);
        end
        vvval(1:Nvu)=1.;
        vvval(vidx(1:vNGLITCH))=vval(1:vNGLITCH);
        filename1=[dirSEED,'/XB.ELYSE.20.BHV.',nameadd];
        filename2=[dirSEED,'/XB.ELYSE.21.BHV.',nameadd];
        filename3=[dirSEED,'/XB.ELYSE.22.BHV.',nameadd];
        filename4=[dirSEED,'/XB.ELYSE.23.BHV.',nameadd];
        filename5=[dirSEED,'/XB.ELYSE.24.BHV.',nameadd];
        mkmseed_IPGP(filename1,glitchv,vt,2,'onefile');
        mkmseed_IPGP(filename2,spredictor(:,2),vt,2,'onefile');
        mkmseed_IPGP(filename3,vdg',vt,2,'onefile');
        mkmseed_IPGP(filename4,v',zt,2,'onefile');
        mkmseed_IPGP(filename5,vvval',zt,2,'onefile');

        %%%%%% save in mseed W in Du
        Nwu=length(w);
        for i=1:Nwu
            glitchw(i)=w(i)-wdg(i)-spredictor(i,3);
        end
        wwval(1:Nwu)=1.;
        wwval(widx(1:wNGLITCH))=wval(1:wNGLITCH);
        filename1=[dirSEED,'/XB.ELYSE.20.BHW.',nameadd];
        filename2=[dirSEED,'/XB.ELYSE.21.BHW.',nameadd];
        filename3=[dirSEED,'/XB.ELYSE.22.BHW.',nameadd];
        filename4=[dirSEED,'/XB.ELYSE.23.BHW.',nameadd];
        filename5=[dirSEED,'/XB.ELYSE.24.BHW.',nameadd];
        mkmseed_IPGP(filename1,glitchw,zt,2,'onefile');
        mkmseed_IPGP(filename2,spredictor(:,3),zt,2,'onefile');
        mkmseed_IPGP(filename3,wdg,zt,2,'onefile');
        mkmseed_IPGP(filename4,w,zt,2,'onefile');
        mkmseed_IPGP(filename5,wwval,zt,2,'onefile');

        filename3=[dirSEED,'/XB.ELYSE.22.BHZ.',nameadd];
        filename4=[dirSEED,'/XB.ELYSE.23.BHZ.',nameadd];
        mkmseed_IPGP(filename3,zzdg,zt,2,'onefile');
        mkmseed_IPGP(filename4,zz,zt,2,'onefile');

        filename3=[dirSEED,'/XB.ELYSE.22.BHN.',nameadd];
        filename4=[dirSEED,'/XB.ELYSE.23.BHN.',nameadd];
        mkmseed_IPGP(filename3,xxdg,zt,2,'onefile');
        mkmseed_IPGP(filename4,xx,zt,2,'onefile');

        filename3=[dirSEED,'/XB.ELYSE.22.BHE.',nameadd];
        filename4=[dirSEED,'/XB.ELYSE.23.BHE.',nameadd];
        mkmseed_IPGP(filename3,yydg,zt,2,'onefile');
        mkmseed_IPGP(filename4,yy,zt,2,'onefile');

    end % end iseed = 1

    %end
end % end if max(NGLITCH_1)>0
end
end
end
end
end 
end % end loop over sols

