% version du 26 octobre 2020
% révision complète par Grégory Sainton.
% simplification massive du code sur des cas très simples

%% DEPENDANCE : 
% DATA:
%    données VBB à 2Hz : vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg
%    réponse instrum   : ./metadata/*
%    FIR filters             : ./FIR/*

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
close all
clc

%% Input parameters

solstart = 82; 
solend  = 82;

%% Read input the data 
% Path to Dropbox
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

sensor = "VBB"; 
mode   = "VEL";

%% plot options
plot_data = 1;
plot_green = 1;

verbose = 1

%% Entry signal parameters

fc=0.005;
sfc=num2str(fc);

dt=.5;
%  attention, dans certains cas ( par exemple iday = 421), 
%  il y a des decallage de plus de 10 msec a expliquer.
%  on reduit ici la tolerance à 10% du samplin
% dt_tol=10e-3;
dt_tol=dt/10;


%% Glitch windows parameters

% Size of the final Green functions
ND2=15*floor(0.5/dt);
NDsum=6*ND2;
NDS=9*ND2;



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


%% Dip and Azimuth + rotation matrix

% Dip et Azi VBB
thetaVBBu=-29.28; phiVBBu=135.11;
thetaVBBv=-29.33; phiVBBv=15.04;
thetaVBBw=-29.61; phiVBBw=254.96; 

% Dip et Azi SP
thetaSPu = -89.9 ; phiSPu = 285.0 ;
thetaSPv = 0.0      ; phiSPv  = 105.2;
thetaSPw = 0.0     ; phiSPw = 345.3;

if sensor == "VBB"
    recompVBB = rotate2zne_mat(thetaVBBu, phiVBBu, thetaVBBv, phiVBBv, thetaVBBw, phiVBBw);
    recomp = recompVBB;
else
    recompSP = rotate2zne_mat(thetaSPu, phiSPu, thetaSPv, phiSPv, thetaSPw, phiSPw);
    recomp = recompSP;
end

for iev=solstart:solend
                                                                                          % condition iextract == 1 supprimée
    igo = 1;                                                                          % set parameters to selected sols
    iday=iev;
        sday=num2str(iday);
    extr='_sol_';
    extr0=strcat('sol',sday);
    
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
    
    close all
    %clearvars -except sensor dirU dirV dirW dirALL dirDROP iday fc sfc recomp epsi dt dt_tol plot_data...
    %  icorrect irot Geom Nhar ifilter ideriv istep valthreshold Tstart Tend extr iday iextract ...
    %    type isoe iyrs imos idas ihrs imns ises iyrf imof idaf ihrf imnf isef ilee iev extr0 dirSEED Nevent igo iseed iequal

    %% lecture des poles et zeros
    
    if sensor == "VBB"
    
        [ugain, uzero, upole, units, uA0]=read_resp_v2('RESP.XB.ELYSE.02.BHU_new');
        [vgain, vzero, vpole, units, vA0]=read_resp_v2('RESP.XB.ELYSE.02.BHV_new');
        [wgain, wzero, wpole, units, wA0]=read_resp_v2('RESP.XB.ELYSE.02.BHW_new');
        C2V_VBB = 335544.0 ;  % from dataless
    elseif sensor == "SP"
        [ugain uzero upole units uA0]=read_resp_v2('RESP.XB.ELYSE.65.EHU');
        [vgain vzero vpole units vA0]=read_resp_v2('RESP.XB.ELYSE.65.EHV');
        [wgain wzero wpole units wA0]=read_resp_v2('RESP.XB.ELYSE.65.EHW');
        C2V_SP = 671088.0;     % from dataless
    end
    Uz=uzero;
    Up=upole;
    Vz=vzero;
    Vp=vpole;
    Wz=wzero;
    Wp=wpole;
      
    %% Read the data
    sday=num2str(iday);
    srot=num2str(irot);
    sthres=num2str(valthreshold);
    sstep=num2str(istep);
    sepsi=num2str(epsi);
    
    filedatau=[dirU,'/sol',sday,'.mat'];
    filedatav=[dirV,'/sol',sday,'.mat'];
    filedataw=[dirW,'/sol',sday,'.mat'];
    

    if exist(filedatau) ~= 0 & exist(filedatav) ~= 0 & exist(filedataw) ~= 0
        % U
        load(filedatau);
        t1=t; ts1=ts; tlmst1=tlmst;
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
    [tz iz]=min(Tfin)
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

        %% Calibration of time series

        %% Remove the Fourier part
        % Fill the gap with zero after Fourier removal

        %%%% u component
        [TTu,bbTTu]=FourierleastsquarewithNaN(utlmst, u, Nhar);
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

        %% equalize the gain ( Here NaN must be removed)
        if iequal == 1
            equalize
        end
        % Suppression de la condition de dérivation

        %% Take only the event interval
        if iextract == 1                           

            [a iau]=min(abs(ut-Tstart));
            [b ibu]=min(abs(ut-Tend));
            [a iav]=min(abs(vt-Tstart));
            [b ibv]=min(abs(vt-Tend));
            [a iaw]=min(abs(wt-Tstart));
            [b ibw]=min(abs(wt-Tend));
            [ut,utlmst,u]=extractanderase0(iau,ibu,ut,utlmst,u);
            [vt,vtlmst,v]=extractanderase0(iav,ibv,vt,vtlmst,v);
            [wt,wtlmst,w]=extractanderase0(iaw,ibw,wt,wtlmst,w);
        else

            Tstart=min([ut(1),vt(1),wt(1)]);
            Tend=max([ut(end),vt(end),wt(end)]);
        end

        %% generate Z axis
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

        % Ajout de la figure en z sur le figure 60 tracée dans la fonction equalize
        subplot(1,4,4)
        plot(ut, uuz, 'k', "DisplayName", "z component")
        legend
        grid on

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
        vdg=v;     % deglitch5 et deglitch7
        wdg=w;
        uuzdg=uuz;

        %% Green function generation
        % Call get_green_func function to create synthetic functions
        % See inside the function for more details
        green_fn_list = get_green_func(sensor, mode, resp, PFO5, PFO2, dt, ND2, NDS, plot_green, verbose);
        
        
        %% Localization of the maxima
        
        vec_data = [u, v, w];
        vec_time = [ut, vt, wt];
        vec_label = ['u', 'v', 'w'];
        vec_struct = []
        
        struct_data = struct;
        
        %loop over axes
        for i_ax=1:size(vec_data, 2) 
            struct_data(i_ax).data = vec_data(:, i_ax); 
            struct_data(i_ax).time = vec_time(:, i_ax); 
            struct_data(i_ax).label = vec_label(:, i_ax); 
            
            prominence = std(vec_data(:, i_ax));                                                                      

            TF_localmin = islocalmin(vec_data(:, i_ax), 'SamplePoints', vec_time(:, i_ax), "MinProminence", prominence);
            TF_localmax = islocalmax(vec_data(:, i_ax), 'SamplePoints', vec_time(:, i_ax), "MinProminence", prominence);

            NblocalMax = find(TF_localmax == 1);
            NblocalMin = find(TF_localmin == 1);
            disp(["Axis: " vec_label(i_ax)])
            disp(["Nb local Min: " length(NblocalMin) " with prominence" prominence])
            disp(["Nb local Max: " length(NblocalMax) " with prominence" prominence])
            
            struct_data(i_ax).peaks = union(NblocalMax , NblocalMin);
            
            
            %ut_d = datetime(vec_time(:, i_ax), 'ConvertFrom', 'datenum');
            ut_d = vec_time(:, i_ax);
            tmp_vec = vec_data(:, i_ax);
            if plot_data ==1
                figure(100+10*i_ax)
                subplot (3,1,1)
                plot(ut_d,vec_data(:, i_ax), 'DisplayName', vec_label(i_ax))
                hold on
                plot(ut_d(TF_localmin), tmp_vec(TF_localmin),'gO',  'DisplayName', "local minima")
                title("Look for local minima")
                grid
                legend

                subplot (3,1,2)
                plot(ut_d,vec_data(:, i_ax), 'DisplayName', vec_label(i_ax))
                hold on
                plot(ut_d(TF_localmax), tmp_vec(TF_localmax),'b*',  'DisplayName', "local maxima")
                title("Look for local maxima")
                grid
                legend

                subplot (3,1,3)
                tmp_vec = abs(vec_data(:, i_ax));
                plot(ut_d, tmp_vec, 'k', 'DisplayName', strcat("abs(",vec_label(i_ax), ")"))
                hold on
                plot(ut_d(TF_localmin),tmp_vec(TF_localmin),'gO',  'DisplayName', "local minima")
                plot(ut_d(TF_localmax),tmp_vec(TF_localmax),'b*',  'DisplayName', "local maxima"')
                title("Absolute value of the signal with the local maxima")
                grid
                legend
            end
        end
        clear tmp_vec
        
        %% Plot any glitch
        figure(200)
        hold on
        grid on
        i_ax  = 1;
        for i=5:1000
            glitch_num = i;
           
            index = struct_data(i_ax).peaks(glitch_num);
            plot(struct_data(i_ax).data(index-ND2:index+NDS))
        end   
        xlabel("VEL glitch")
        title("Overplot glitches")
        
      
        
    end
end
end
end
end 
end % end loop over sols

%% Definition of the function get_green_func
function green_fn_list = get_green_func(sensor, mode, resp, PFO5, PFO2, dt, ND2, NDS, plot_green, verbose)
    % Function which estimate the green functions of a step function
    % 
    % We are assuming that a glitch can be modeled as a step in
    % acceleration. 
    % 1> Creation of a transfert function (TF )using zpk matlab function : 
    %                https://www.mathworks.com/help/control/ref/zpk.html
    % 2> Apply this TF to a step function
    %
    % 3> Depending on the content of type_Green the zero list is modified
    % ----------
    % Parameters:
    %       @sensor       : VBB or SP
    %       @mode        : "VEL" or "POS"
    %       @resp          : array with the path of the dataless (contaning the transfert functions)
    %       @PFO5         : FIR coeff with gain 5
    %       @PFO2         : FIR coef with gain 2
    %       @dt              : sampling rate of the signal 
    %       @ND2           : # of point before max of the GF
    %       @NDS           : # of point after the max of the GF     
    %       @plot_green  : 0 or 1 (1 if you want to plot the green functions)
    %       @verbose      : 0 or 1 (1 if you want informations about transfer functions and so on... )
    % -----------
    % Outputs
    %       @green_fn_list: array of array where each array is a green function 
    %                              according to the list type_green
    %
    % ------------
    % Dependance
    %       decimeFPGA : decimation using the FIR coefs filters
    
    % Parameters 
    type_green = ["G", "Gprime", "P", "Pprime"];   
    
    NS = 1024;
    dti=dt/10;        % Output sampling rate
    NSi=10*NS;   
     
    % Read of the transfert functions
    %        TF are saved in the directory "metadata"
   
     if sensor == "VBB" && mode == "VEL"
            [ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.vbb_vel_hg);
     elseif sensor == "VBB" && mode == "POS"
          [ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.vbb_pos_hg);
     elseif sensor == "SP"
          [ugain, uzero_init, upole_init, units, uA0]=read_resp_v2(resp.sp_hg);
     end
      
    if verbose == 1
        disp(strcat("For the sensor ", sensor, "with mode", mode, ": "))
        disp(["ugain= " ugain])
        disp(["uA0" uA0])
        disp(["uzero= " uzero_init]) 
        disp(["upole= " upole_init])
        disp(["Units = " units])
    end
     
     green_fn_list = [];                                 % Output var with the list of the green functions
     time_20sps=[0:NSi-1]*dti;                    % Time vector
    
     % Modify pole list if VEl or POS
     if mode == "VEL"
         upole = [0 upole_init];
     else 
         upole = upole_init;
     end
   
   % Estimate the first Green function  
   %
   sys = zpk(uzero_init, upole, ugain*uA0);                              % Ca ne marche pas en introduisant le Gain dans zpk
   %sys = zpk(uzero, upole, ugain, dti);                                  % Les gains sont incohérents                            
   sys.TimeUnit = "seconds";
   
    if plot_green == 1 
       figure(300)
       bode(sys)
       h = gcr;
       setoptions(h,'FreqUnits','Hz');
       grid on 
   end
  
   % Apply a step function to the TF
   green_20sps= step(sys, time_20sps);   
   
   % Center the function 
   NF=floor(NSi/2+1);                                                           % middle point of the vector
   green_20sps_shift = [zeros(1, NF), green_20sps'];              % Shift the vector to center the max of the Green func
   green_20sps_shift = green_20sps_shift(1:NSi);
   
   % Append in a list
   green_fn_list_20 = [green_20sps_shift];
   
   % Calculation of the derivates of G
   % To avoid causality issues, calculation is made using the TF derivation
   % f^{(k)}_m = (imw)^kf_m$ avec $f_m = \frac{1}{T}\int_0^T e^{imwt}f(t)dt  
   for deriv=1:3
       omega(1:NF)=[0:NF-1]*1i/NS/dti*2*pi;
       tf=fft(green_20sps_shift);  
       tf(1:NF)=tf(1:NF).*(omega(1:NF));
       green_20sps_shift=ifft(tf,'symmetric');
       green_fn_list_20 = [green_fn_list_20; green_20sps_shift];
   end
   
  [numRows,numCols] = size(green_fn_list_20);
  
   if plot_green == 1 
         figure(301);
   end 
   for i=1:numRows
       
        % decimation from 20sps to 2sps
        %    TODO : add parametrization to avoid such hard coding
        green_4sps  =  decimeFPGA(green_fn_list_20(i,:), 5, PFO5, sum(PFO5));
        green_2sps  =  decimeFPGA(green_4sps, 2, PFO2, sum(PFO2));

        [amp_max_green_2sps, idx_max_2sps] = max(abs(green_2sps));
        if i==1
            idx_ref = idx_max_2sps;
        end
        green = green_2sps(idx_ref-ND2: idx_ref+NDS);  % normalisation avec le maximum
                                                                                           % On garde 15 points avant et 135 après
        % Append the function to the output list                                                                                      
        green_fn_list = [green_fn_list ; green];

        if plot_green == 1
            cmap = hsv(numRows); 
            sgtitle("Synthetic functions")
            subplot(numRows, 1, i )
            plot(green, "g", "DisplayName" , type_green(i), 'Color', cmap(i,:))
            ylabel("D.U./m.s^{-2}")
            xlabel("Time in seconds")
            xlim([0 150])
            legend
            grid on
            title(type_green(i))
            
        end            % end plot option 
   end                 % loop over function 
   clear PFO2 PFO5
end                    % end get_green_func 
