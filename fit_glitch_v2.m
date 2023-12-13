%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% name     : fit_glitch.m
% author   : initial algorithm from Philippe Lognonné
% revision : Grégory Sainton  
% 
% purpose: Detect and remove glitches on seismic signals from 
%          SEIS seismometers
%     
%
% For the moment, the code does not take into account partial sol of data
% Input data are at 20 Hz sps and output data are at 2 sps
% 
% C'est le programme principal.
% Toute la partie préprocessing (ie ligne 458 environ), c'est le même
% processing que celui proposé par Philippe dans son code. La partie
% spécifique commence après.
% 
% La partie "Fit parameters" est celle où on règle les paramètres de fit des glitches
% La partie "Plot and verbose options" -> Tout est dans le titre, ce sont des flags


%% DEPENDANCES : 
% DATA:
%    VBB Data @ 2Hz           : vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg
%    Instrumental responses   : ./metadata/*
%    FIR filters              : ./FIR/*
%
% FUNCTIONS : 
%   All the functions are saved in subroutines_PL 
%                               or subroutines_GS
%

clear all
close all
clc
warning('off','all')

%% Input parameters

solstart = 999; 
solend  = 1000;

%% Read input the data 
% Path to Dropbox
   
dirDATA = '/Users/Greg/Dropbox (IPGP)/';
dirDROP = '/Users/Greg/CODE/';
dirCAT =  'VBB_Glitches_analysis/'
%dirDROP         = '/Users/gr3g/Dropbox (IPGP)/';
dirCodeDeglitch = strcat(dirDROP,"deglitch/");
outputdir       = strcat(dirDATA, dirCAT);

% Path to 2 sps VEL High Gain VBB data for the 3 components, U, V, W
%origindir = 'vbb_data_for_philippe';
origindir = 'SEIS_ON_MARS_IPGP/DAILY_INSIGHT_DATA/';
dirU            = [dirDATA origindir 'matlab_data/vel_sci/2sps/hg/vbbue'];  % Read only dir
dirV            = [dirDATA origindir 'matlab_data/vel_sci/2sps/hg/vbbve'];  % Read only dir
dirW            = [dirDATA origindir 'matlab_data/vel_sci/2sps/hg/vbbwe']; % Read only dir

% Dir to the FIR filter files
dir_FIR         = [dirDROP  'deglitch/FIR/'];
PFO5            = load([dir_FIR 'PFO_div5.txt']);      % FIR coefs to decimate by 5
PFO4            = load([dir_FIR 'PFO_div4.txt']);      % FIR coefs to decimate by 5
PFO2            = load([dir_FIR 'PFO_div2.txt']);      % FIR coefs to decimate by 2
PFO_list        = {PFO2, PFO4, PFO5};

% Dir to the metadata
dir_resp        = [dirDROP  'deglitch/metadata/'];
resp.vbb_vel_hg = [dir_resp  'RESP.XB.ELYSE.02.BHU_new'];
resp.vbb_pos_hg = [dir_resp  'RESP.XB.ELYSE.00.LMU'];
resp.sp_hg      = [dir_resp  'RESP.XB.ELYSE.65.EHU'];

%% Primary set of parameters

iseed=0;
icorrect=0;
ifilter=0;
irot=1;

Nhar=12;                                  % Number of harmonics to detrend the signal (related to "FourierleastsquarewithNaN")
iextract=0;

iequal=1;

sensor = "VBB"; 
mode   = "VEL";


%% -------   Fit parameters -------
%
 
typefit = "gp";                          % "g" for Green or "gp" for Green & Precursor
lagrangemult = 1;                    % add Lagrange multipliers to constrain the first and the list point of the fit 
polyfit = 0  ;                             % if =1 -> Entangled glitches are fitted with a single function otherwise, 
                                                %        the glitches into a cluster are fitted one by one in descending   amplitude

prom_factor = 1.5;                   % prominence factor used to detect the local extrema in the function find_local_maxima
cut_fact = 0.6;                          

param_fit_glitch = 1;                % parameter to allow the glitch fitting or not

dt=.5;                                       % Frequency of the signal
dist_min_glitch = 1.25/dt;         % threshold to merge extrema from different axes into a single one
radius = 20 ;                             % minimal distance to consider an extremum single
radius_min =0 ;
radius_max = radius;

fc=0.005;
sfc=num2str(fc);                                                                 
dt_tol=dt/10;

% Size of the final Green functions (GF)
ND2=15*floor(0.5/dt);                % Number of points before the max amplitude of the GF 
NDS=floor(4*ND2);                     % Number of points after the max (Philippe's value is %NDsum=6*ND2;)
                                                                 

% Cut to eliminate "rebound" of the glitches considered as
% extrema by the functions islocalmin and islocalmax
epsG = 4;
gapG = 20;

% Final cut in variance
var_thres = 0.04;   % Variance threshold to cut high variance glitches
%var_thres = 0.11;   

%% plot and verbose options
plot_data         = 0;                        % figure(60)
plot_green       = 0;                        % figure(300) and figure (301)
plot_extrema   = 0;                        % figure(101) to figure(103) 
plot_fit             = 0;                        % figure(400) to figure(420)
verbose            = 1;                       % Display explanation during the process
plot_fit_clusters  = 0;                      % To plot all the glitches on by one -> VERY TIME CONSUMING
save_ts_data     = 1;
save_param      = 1;
%% Output directories with deglichted files in both MSEED and MATformat
if icorrect == 0
    dirALL  =[dirDATA 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGFP0_var1'];
    dirSEED=[dirDATA 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGSP0'];
else
    if istep == 1     %one step of deglitching
        dirALL=[dirDATA 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGF1'];
    elseif istep == 2 %two steps of deglitching
        dirALL=[dirDATA 'DeglitchedData_PL/vel_sci/2sps/hg/vbbDGF2'];      
    end
end

%% Dip and Azimuth + rotation matrix

% Dip et Azi VBB  -> Should be taken from dataless
thetaVBBu=-29.28; phiVBBu=135.11;
thetaVBBv=-29.33; phiVBBv=15.04;
thetaVBBw=-29.61; phiVBBw=254.96; 

% Dip et Azi SP -> Should be taken from dataless
thetaSPu = -89.9  ; phiSPu = 285.0 ;
thetaSPv = 0.0      ; phiSPv  = 105.2 ;
thetaSPw = 0.0     ; phiSPw = 345.3 ;

if sensor == "VBB"
    recompVBB = rotate2zne_mat(thetaVBBu, phiVBBu, thetaVBBv, phiVBBv, thetaVBBw, phiVBBw);
    recomp    = recompVBB;
else
    recompSP  = rotate2zne_mat(thetaSPu, phiSPu, thetaSPv, phiSPv, thetaSPw, phiSPw);
    recomp    = recompSP;
end

%% poles and zeros reading

    if sensor == "VBB" && mode == "VEL"

        [ugain, uzero, upole, unitsu, uA0]=read_resp_v2('RESP.XB.ELYSE.02.BHU_new');
        [vgain, vzero, vpole, unitsv, vA0]=read_resp_v2('RESP.XB.ELYSE.02.BHV_new');
        [wgain, wzero, wpole, unitsw, wA0]=read_resp_v2('RESP.XB.ELYSE.02.BHW_new');
        C2V_VBB = 335544.0 ;  % from dataless 
    elseif sensor == "VBB" && mode == "POS"     
       [ugain, uzero, upole, unitsu, uA0]=read_resp_v2('RESP.XB.ELYSE.00.LMU');
       [vgain, vzero, vpole, unitsv, vA0]=read_resp_v2('RESP.XB.ELYSE.00.LMV');
       [wgain, wzero, wpole, unitsw, wA0]=read_resp_v2('RESP.XB.ELYSE.00.LMW');
    elseif sensor == "SP"
        [ugain, uzero, upole, unitsu, uA0]=read_resp_v2('RESP.XB.ELYSE.65.EHU');
        [vgain, vzero, vpole, unitsv, vA0]=read_resp_v2('RESP.XB.ELYSE.65.EHV');
        [wgain, wzero, wpole, unitsw, wA0]=read_resp_v2('RESP.XB.ELYSE.65.EHW');
        C2V_SP = 671088.0;     % from dataless
    end
    
    
%% LOOP OVER SOLS

for iev=solstart:solend
    disp(["Sol number; ", iev])                                                       % condition iextract == 1 supprimée
    
    tsfilesave       = strcat(outputdir,"TS_",num2str(iev),".mat");
    
    igo = 1;                                                                          % set parameters to selected sols
    iday=iev;
    sday=num2str(iday);
    extr='_sol_';
    extr0=strcat('sol', sday);
    
    if igo == 1

        if iday  >=74   && iday <= 167  % data were already @ 2sps
            dirU=[dirDATA 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbu'];
            dirV=[dirDATA 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbv'];
            dirW=[dirDATA 'vbb_data_for_philippe/matlab_data/vel_sci/2sps/hg/vbbw'];
        elseif iday >= 168               % data decimated from 20Hz to 2Hz 
            dirU            = [dirDATA origindir '/matlab_data/vel_sci/2sps/hg/vbbue']; % Read only dir
            dirV            = [dirDATA origindir '/matlab_data/vel_sci/2sps/hg/vbbve']; % Read only dir
            dirW            = [dirDATA origindir '/matlab_data/vel_sci/2sps/hg/vbbwe']; % Read only dir
        end
    
        %close all

        Uz=uzero;
        Up=upole;
        Vz=vzero;
        Vp=vpole;
        Wz=wzero;
        Wp=wpole;

        %% Read the data
        sday=num2str(iday);
        srot=num2str(irot);
        %sthres=num2str(valthreshold);
    
        filedatau=[dirU,'/sol',sday,'.mat']
        filedatav=[dirV,'/sol',sday,'.mat']
        filedataw=[dirW,'/sol',sday,'.mat']

        if exist(filedatau) ~= 0 && exist(filedatav) ~= 0 && exist(filedataw) ~= 0
            % U
            load(filedatau);
            t1=t; ts1=ts; tlmst1=tlmst;
            % V
            load(filedatav);
            t2=t; ts2=ts; tlmst2=tlmst;
            % W
            load(filedataw);
            t3=t; ts3=ts; tlmst3=tlmst;

            if ~isempty(ts1) && ~isempty(ts2) && ~isempty(ts3) 

            % check missing data and synchronize on common start of u,v,w

                kk1=1;
                while(ts1(kk1) == 0 || ts1(kk1) ~= ts1(kk1)) 
                    kk1=kk1+1;
                end
                kk2=1;
                while(ts2(kk2) == 0 || ts2(kk2) ~= ts2(kk2)) 
                    kk2=kk2+1;
                end
                kk3=1;
                while(ts3(kk3) == 0 || ts3(kk3) ~= ts3(kk3))
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
                    % detect any data interuption on u,v,w
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
                    [tz, iz]=min(Tfin);
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

                    utlmst=tlmst1(kk1:Nu);
                    ut=t1(kk1:Nu);
                    u=ts1(kk1:Nu);
                    vtlmst=tlmst2(kk2:Nv);
                    vt=t2(kk2:Nv);
                    v=ts2(kk2:Nv);
                    wtlmst=tlmst3(kk3:Nw);
                    wt=t3(kk3:Nw);
                    w=ts3(kk3:Nw);

                    % force interpolation on the longest components in case of holes
                    Nu=length(ut);
                    Nv=length(vt);
                    Nw=length(wt);
                    Ncase=[Nu, Nv, Nw];
                    [Nmax, ic]=max(Ncase);
                    Ncase=Ncase-Nmax;
                    if max(abs(Ncase)) > 0
                        warning1='ATTENTION LES COMPOSANTES U,V,W N ONT PAS LE MEME NOMBRE DE POINTS';
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
                    % flag start and end
                    urawstart=u(1);urawend=u(end);
                    vrawstart=v(1);vrawend=v(end);
                    wrawstart=w(1);wrawend=w(end);

                    %% Calibration of time series %%

                    %% 1. Remove the Fourier part
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
                    
                    
                    %% Keep unequalized series (used in the final part to saved the data)
                    u_orig = u;
                    v_orig  = v;
                    w_orig = w;
                    
                    %% 2. Equalize the gain ( Here NaN must be removed)
                    if iequal == 1
                        equalize
                    end
                    % Suppression de la condition de dérivation

                    %% Take only the event interval
                    if iextract == 1                           

                        [~, iau]     = min(abs(ut-Tstart));
                        [~, ibu]     = min(abs(ut-Tend));
                        [~, iav]     = min(abs(vt-Tstart));
                        [~, ibv]     = min(abs(vt-Tend));
                        [a, iaw]     = min(abs(wt-Tstart));
                        [b, ibw]     = min(abs(wt-Tend));
                        [ut, utlmst, u]= extractanderase0(iau,ibu,ut,utlmst,u);
                        [vt, vtlmst, v]= extractanderase0(iav,ibv,vt,vtlmst,v);
                        [wt, wtlmst, w]= extractanderase0(iaw,ibw,wt,wtlmst,w);
                    else

                        Tstart=min([ut(1),vt(1),wt(1)]);
                        Tend=max([ut(end),vt(end),wt(end)]);
                    end

                    %% 3. Generate Z axis
                    uuz=recomp(1,1)*u+recomp(1,2)*v+recomp(1,3)*w;
                    uux=recomp(2,1)*u+recomp(2,2)*v+recomp(2,3)*w;
                    uuy=recomp(3,1)*u+recomp(3,2)*v+recomp(3,3)*w;
                    TTuuz=recomp(1,1)*TTu+recomp(1,2)*TTv+recomp(1,3)*TTw;
                    TTuux=recomp(2,1)*TTu+recomp(2,2)*TTv+recomp(2,3)*TTw;
                    TTuuy=recomp(3,1)*TTu+recomp(3,2)*TTv+recomp(3,3)*TTw;
                    bbTTuuz=recomp(1,1)*bbTTu+recomp(1,2)*bbTTv+recomp(1,3)*bbTTw;
                    bbTTuux=recomp(2,1)*bbTTu+recomp(2,2)*bbTTv+recomp(2,3)*bbTTw;
                    bbTTuuy=recomp(3,1)*bbTTu+recomp(3,2)*bbTTv+recomp(3,3)*bbTTw;

                    % Ajout de la figure en z sur le figure 60 tracée dans la fonction equalize
                    if plot_data ==1
                        subplot(4,1,4)
                        plot(ut, uuz, 'k', "DisplayName", "z component")
                        legend
                        grid on
                    end

                    %%%% End of preprocessing %%%%
       
                    %% Green function generation
                    % Call get_green_func function to create synthetic functions
                    % See inside the function for more details
                    green_fn_list = get_green_func(sensor, mode, resp, PFO_list, dt, ND2, NDS, plot_green, verbose);

                    % Creation of a vector of Weight to overweight the max of the
                    % function and the signal during fit
                    [agdelay, gdelay]=max(abs(green_fn_list(1, :)));

                    dtweight=4.;
                    ntweight=floor(dtweight/dt);
                    weight(1:length(green_fn_list(1,:)))=1.;
                    %weight(1:gdelay/2.)=50.;
                    %weight(gdelay/2.0:gdelay+3.5*ntweight)=100.;
                    weight(1:gdelay+3.5*ntweight)=100.;
                    %weight = abs(green_fn_list(1,:))/(sum(abs(green_fn_list(1,:))));
                    
                    
                    if plot_green ==1
                        figure(302)
                        plot(weight)
                        title("Weight vector")
                        xlabel("Time")
                        ylabel("Weight")
                        grid on 
                        legend
                    end

                    %% Prepare structure of data -> used to save all the data
                     struct_data = struct;

                     vec_data = [u, v, w];
                     vec_time = [ut, vt, wt];
                     vec_lmst = [utlmst, vtlmst, wtlmst];
                     vec_label = ['u', 'v', 'w'];

                     for i_ax=1:3
                        struct_data(i_ax).data = vec_data(:, i_ax); 
                        struct_data(i_ax).time = vec_time(:, i_ax); 
                        struct_data(i_ax).lmst = vec_lmst(:, i_ax); 
                        struct_data(i_ax).label = vec_label(:, i_ax);
                     end
                    clear vec_data vec_time vec_label 
                    
                    %% Localization of the maxima
                    %  -------------------------------
                    % Find maxima
                   
                    struct_data              = find_local_maxima2(struct_data, prom_factor, cut_fact, plot_extrema, verbose);
                    
                    % Merge maxima for each axes
                    full_glitch_vec          = union(struct_data(1).peaks, union(struct_data(2).peaks, struct_data(3).peaks));  
                    max_glitch_per_axis = max([struct_data(:).nbpeaks]);
                    
                    % Merge the indexes to close which are probably the same glitches detected
                    % on the different axes.
                    reduced_glitch_vec  = merge_close_glitches(full_glitch_vec, dist_min_glitch);
                    
                    if verbose == 1
                        disp(["Longueur de reduced_glitch_vec avant cut: ", length(reduced_glitch_vec)])
                    end

                    % Cut to eliminate the maxima considered as glitches by
                    % find_local_maxima2 function and which are actually 
                    % the rebound part of the Green function due to the
                    % transfert function
                   
                     for i=1:length(reduced_glitch_vec)-1  
                         if (abs(reduced_glitch_vec(i)-reduced_glitch_vec(i+1)) >= gapG - epsG) && ...                                % select extrema inside range gapG +/- epsG 
                            (abs(reduced_glitch_vec(i)-reduced_glitch_vec(i+1))< gapG + epsG)                   
                             if struct_data(i_ax).data(reduced_glitch_vec(i))*struct_data(i_ax).data(reduced_glitch_vec(i+1))<0  % rebound must have opposite sign with 
                                                                                                                                                                           % primary signal
                                 reduced_glitch_vec(i+1)=NaN;                                                                                              % set to NaN because convienient to remove
                             end
                         end
                     end
                     reduced_glitch_vec = reduced_glitch_vec(~isnan(reduced_glitch_vec));                                        % Finally remove the NaN created above
                     
                     
                    if verbose == 1
                        disp(["Longueur de reduced_glitch_vec après cut: ", length(reduced_glitch_vec)])
                        disp(strcat("Merge the extrema with gap less than ", num2str(dist_min_glitch)))
                        disp(strcat(num2str(length(full_glitch_vec)), " before"))
                        disp(strcat(num2str(length(reduced_glitch_vec)), " after"))
                    end
                    
                    if length(reduced_glitch_vec)>=1
                    
                        % Create a list of mono-glitches -> they are easy to
                        % fit because no entanglement with another glitches
                        single_extrema_list = get_single_extrema(reduced_glitch_vec, radius);
                        results_fit         = zeros(length(single_extrema_list), 1, 3);

                        if verbose == 1
                            disp(strcat("There are ", num2str(length(single_extrema_list)), ...
                                 " single extrema to fit alone."))
                        end

                        % From the list of single glitches one can extract the
                        % list of polyglitches.
                        mask                                  = ismember(reduced_glitch_vec, single_extrema_list);
                        multi_extrema_list                = reduced_glitch_vec;
                        multi_extrema_list(mask)      = [];
                        clear mask

                        if verbose == 1
                            disp(strcat("There are  ", num2str(length(multi_extrema_list)), ...
                                 " multiple extrema to fit in clusters."))
                         end

                        %% Fitting part of the glitches
                        if param_fit_glitch == 1                                    % Flag to skip the fit

                            % FIRST, We are dealing with the single glitches
                            if length(single_extrema_list)>=1

                                index_list_2_deglitch         = single_extrema_list;
                                [struct_data, Glitch_cat]     =  fit_list_of_single_gliches_v2(struct_data, index_list_2_deglitch, ...
                                                                                                      max_glitch_per_axis, green_fn_list, ...
                                                                                                      typefit, lagrangemult, weight, var_thres, ...
                                                                                                      ND2, NDS, verbose);
                                
                                %Plot of the synthetic glithes
%                                     figure(2000+iev)
%                                     axes = [];
%                                     for i_ax=1:3
%                                         axes(i_ax)= subplot(3,1, i_ax);
%                                         plot(struct_data(i_ax).data, "k" ,'DisplayName', strcat("Original signal ", struct_data(i_ax).label))
% 
%                                         hold on 
%                                         plot(struct_data(i_ax).synth_glitch,"r" ,'DisplayName', strcat("Synthethic glitch on ", struct_data(i_ax).label))
%                                         plot(struct_data(i_ax).data_dg,     "g", 'DisplayName', strcat("Deglitched signal on ", struct_data(i_ax).label))
%                                         xlabel("Index")
%                                         ylabel("Amplitude")
%                                         grid on
%                                         
%                                          %for i=1:length([substruct.index])
%                                          for i=1:length([Glitch_cat.index])
%                                              scatter(Glitch_cat(i).index, Glitch_cat(i).(strcat("maxAmpFit_",struct_data(i_ax).label)) ,"c*")
%                                          end
%                                          clear i
% 
%                                         clear idx clust_idx gl
%                                         legend('Original signal', 'Glitches signal', 'Deglitched signal')
%                                         legend 
%                                     end % loop over axes
%                                     linkaxes(axes, 'xy');
% 
%                                     clearvars i_ax
                            end

                            % SECOND, We are dealing with the multiple glitches
                             if length(multi_extrema_list)>=2
                                clusters_of_glitches                          = clusterize_extrema(multi_extrema_list, radius_min, radius_max);
                            
                                if polyfit == 0
                                    % Glitches are of a given cluster are fitted in
                                    % descending order of amplitude.
                                    [struct_data, Glitch_cat]     = fit_list_of_multi_gliches_v4(struct_data, clusters_of_glitches, ...
                                                                                                           max_glitch_per_axis, green_fn_list, ...
                                                                                                            dt, typefit, lagrangemult, ...
                                                                                                           weight, ND2, NDS, Glitch_cat, var_thres, ...
                                                                                                           verbose, plot_fit_clusters);

                                else
                                    % Glitches are of a given cluster are fitted in
                                    % one time.
                                    [struct_data,  Glitch_cat]     = fit_list_of_multi_gliches(struct_data, clusters_of_glitches, ...
                                                                                                       max_glitch_per_axis, green_fn_list, ...
                                                                                                        dt, typefit, lagrangemult, ...
                                                                                                       ND2, NDS, Glitch_cat, verbose);

                                end
                            else
                                disp("No cluster of glitches in this SOL")
                                
                            end
                            % TRIGGER PART TO CLEAN GLITCHES 
                            if length(Glitch_cat)>=1
                                T                           = struct2table(Glitch_cat);
                                sortedTabGlitchCat = sortrows(T, 'index');
                                Glitch_cat               = table2struct(sortedTabGlitchCat);     
                                AllGlitchesFitted      = Glitch_cat;
                                % Select the glitches where the VARIANCE reduction
                                % (before/after fit) is lower than "var_thres"

                                %subGlitch_cat = sortedTabGlitchCat((sortedTabGlitchCat.var_tot_ratio_u <= var_thres) | ...
                                %                                    (sortedTabGlitchCat.var_tot_ratio_v <= var_thres) | ...
                                %                                    (sortedTabGlitchCat.var_tot_ratio_w <= var_thres), :);
                                                                
                                subGlitch_cat = sortedTabGlitchCat((sortedTabGlitchCat.isglitch_u == 1) | ...
                                                                    (sortedTabGlitchCat.isglitch_v == 1) | ...
                                                                    (sortedTabGlitchCat.isglitch_w == 1), :);
                                                        
                               if verbose == 1
                                    disp(strcat("    Variance reduction cut @ ", num2str(var_thres), ...
                                                " It remains ", num2str(height(subGlitch_cat)), "  glitches."))
                                end 
                                substruct = table2struct(subGlitch_cat);
                                Glitch_cat = substruct;

                                % Select the glitches where the AMPLITUDE is larger
                                % than a standard deviation to avoid keeping to
                                % much noise.
                                ampli_cut = 0;
                                if ampli_cut == 1
                                    amp_thres_u = mean(struct_data(1).data) + std(struct_data(1).data);
                                    amp_thres_v = mean(struct_data(2).data) + std(struct_data(2).data);
                                    amp_thres_w = mean(struct_data(3).data) + std(struct_data(3).data);
                                    subGlitch_cat = subGlitch_cat(abs(subGlitch_cat.data_u) > amp_thres_u | ...
                                                                                     abs(subGlitch_cat.data_v) > amp_thres_v | ...
                                                                                     abs(subGlitch_cat.data_w) > amp_thres_w, :);

                                    if verbose == 1
                                        disp(strcat("Amplitude cut @, it remains ", num2str(height(subGlitch_cat)), "  glitches."))
                                    end
                                    substruct = table2struct(subGlitch_cat);
                                    Glitch_cat = substruct;
                                end 
%                                 %% Plot the result of the fits
%                                 %  ------------------------------------------------
% 
                                if plot_fit == 1

                                    %Plot of the synthetic glithes
                                    figure(1000+iev)
                                    axes = [];
                                    for i_ax=1:3
                                        axes(i_ax)= subplot(3,1, i_ax);
                                        plot(struct_data(i_ax).data,"k" ,'DisplayName', strcat("Original signal ", struct_data(i_ax).label))

                                        hold on 
                                        plot(struct_data(i_ax).synth_glitch,"r" ,'DisplayName', strcat("Synthethic glitch on ", struct_data(i_ax).label))
                                        plot(struct_data(i_ax).data_dg,     "g", 'DisplayName', strcat("Deglitched signal on ", struct_data(i_ax).label))
                                        xlabel("Index")
                                        ylabel("Amplitude")
                                        grid on

                                         %for i=1:length([substruct.index])
                                         for i=1:length([Glitch_cat.index])
                                             scatter(Glitch_cat(i).index, struct_data(i_ax).data(Glitch_cat(i).index),"c*")
                                         end
                                         clear i

                                        clear idx clust_idx gl
                                        legend('Original signal', 'Glitches signal', 'Deglitched signal')
                                        legend 
                                    end % loop over axes
                                    linkaxes(axes, 'xy');

                                    clearvars i_ax
                                end
                             end % param_fit_glitch =1
                        end
                    end
                    
                        %% SAVING PART 
                        if save_param == 1
                            % Results are saved into files. 
                            glitchfile = strcat(outputdir, "Catalogue/", "SOL",num2str(iev, "%04d"), "_Glitches_cat.mat");
                            save (glitchfile, 'Glitch_cat')

                            %greenfunc = strcat(dirCodeDeglitch, "Output/", "Green_functions.mat");
                            %save (greenfunc, 'green_fn_list')

                            if save_ts_data == 1
                               save(tsfilesave, 'struct_data')
                            end

                            % extract data from structure for saving part.
                            udg = struct_data(1).data_dg; ugl = struct_data(1).synth_glitch;
                            vdg = struct_data(2).data_dg; vgl = struct_data(2).synth_glitch;
                            wdg =struct_data(3).data_dg; wgl = struct_data(3).synth_glitch;

                            % Remove equalization effect on amplitude of the deglitched data
                            udg=remove_equalizing(udg, Uz, Up, uA0, ugain, Uz, Up, uA0, ugain, dt);
                            vdg=remove_equalizing(vdg, Vz, Vp, vA0, vgain, Uz, Up, uA0, ugain, dt);
                            wdg=remove_equalizing(wdg, Wz, Wp, wA0, wgain, Uz, Up, uA0, ugain, dt);

                            % Remove equalization effect on amplitude of the glitch signal
                            ugl = remove_equalizing(ugl, Uz, Up, uA0, ugain, Uz, Up, uA0, ugain, dt);
                            vgl = remove_equalizing(vgl, Vz, Vp, vA0, vgain, Uz, Up, uA0, ugain, dt);
                            wgl = remove_equalizing(wgl, Wz, Wp, wA0, wgain, Uz, Up, uA0, ugain, dt);

                            outdatafile = strcat(outputdir, "Output/", "SOL",num2str(iev, "%04d"),"_", num2str(var_thres) ,"_outdata.mat");
                            save (outdatafile, 'utlmst', 'ut', 'u', 'udg', 'ugl', ...
                                                    'vtlmst', 'vt', 'v', 'vdg', 'vgl', ...
                                                    'wtlmst', 'wt', 'w', 'wdg', 'wgl')
                            disp(strcat("Glitches catalog of sol ", num2str(iev), " saved at ", outdatafile))

                        end
                end
            else
                disp("One channel is empty")
            end
        else
            disp("Some files are empty")
        end
    else
        disp("NO GO")
    end
end % end loop over sols