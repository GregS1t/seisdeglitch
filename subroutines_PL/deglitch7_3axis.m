%%%%% deglitch data in 3 axis configuration
%%%%% P.Lognonné, 2/27/2019
%%%%% v0
%% Parameters
% INPUT
% ---------
%    ut : timing
%     u  raw data
%     NS Length of the Green function in points
%     dt Sampling time
%     uz zero of Transfer function
%     up pole of Transfer function
%     NZERO  number of glitch for search
%     valthreshold threshold for glitch acceptance in variance reduction
%     dthreshold threshold for glitch acceptance in time offset
% OUTPUT
% -------
%     unew : deglitched signal
%     uglicth: glitch signal
%     Var0:  variance reduction for the input signal
%     Number of Glitch removed
%     val variance recution
%     idx starting index of the glitch
%     Amp inverted parameters ( glitch amplitude, glitch time shift, mean, trend, square drift)
%%%%%%%%% more ( 9/9/2020)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [unew,uglitch,udrift,udiff,u,ustack,ustack0,...
            vnew,vglitch,vdrift,vdiff,v,vstack,vstack0,...
            wnew,wglitch,wdrift,wdiff,w,wstack,wstack0,...
            uuznew,uuzdiff,NGLITCH,Var0,Var1,val,val0,...
            idx,Amp,Delay,Green,ND,im,g0,ND2,NDS,nstack,...
            NG1,NG2,NZEROALL,sprecursor,jdx,AmpSpike,reducA, Var]...
            =deglitch7_3axis(...
                ut,utlmst,udg,u,vt,vtlmst,vdg,v,...
                wt,wtlmst,wdg,w,uuzdg,uuz,recomp,epsi,NS,dt,...
                Uz,Up,Vz,Vp,Wz,Wp,NZERO,valthreshold,dtthreshold,tfthreshold,...
                locm,NLOC, Nderiv, istep, iday, icorrect , igg, iequal)


dirDROP = '/Users/gr3g/Dropbox (IPGP)/';
%dirDROP = '/Users/greg/Dropbox (IPGP)/';
DIRSTACKF = [dirDROP 'DeglitchedData_PL/dirstackF/'];
DIRSTACK20SPS = [dirDROP 'DeglitchedData_PL/dirstack20sps/'];
DIRSTACK_CORRECT_F = [dirDROP 'DeglitchedData_PL/dirstack_correctF/'];

%%%%% remove the very long period to ensure to drift/offset or trend

iplot=1;
ifilter=0;
immax=10;
jdx=[];
AmpSpike=[];
reducA=[];
if ifilter == 2
	[B,A]=butter(4,0.2*dt,'high');         %%%% low pass filter parameter
	udiff(1:length(u))=0.;
	vdiff(1:length(v))=0.;
	wdiff(1:length(w))=0.;
	uuzdiff(1:length(u))=0.;
elseif ifilter == 1
	[B,A]=butter(5,.005*dt,'low');         %%%% low pass filter parameter
	udiff=filtfilt(B,A,udg);                      %%%% remove the trend from the  deglitched input data
	vdiff=filtfilt(B,A,vdg);               %%%% remove the trend from the  deglitched input data
	wdiff=filtfilt(B,A,wdg);               %%%% remove the trend from the  deglitched input data
	uuzdiff=filtfilt(B,A,uuzdg);           %%%% remove the trend from the  deglitched input data
	u=u-udiff;                             %%%% take out the trend 
	v=v-vdiff;                             %%%% take out the trend 
	w=w-wdiff;                             %%%% take out the trend 
	uuz=uuz-uuzdiff;                       %%%% take out the trend 
else
	udiff(1:length(u))=0.;
	vdiff(1:length(v))=0.;
	wdiff(1:length(w))=0.;
	uuzdiff(1:length(u))=0.;
end
s(:,1)=u;
s(:,2)=v;
s(:,3)=w;

%%%%% built Green response from Pole and Zero ( remove the two first null zero)

icase=2;

%% First loop on axes 
for iaxe=1:3

    if iequal == 0
        if iaxe == 1
            up=Up;
            uz=Uz;
        elseif iaxe == 2
            up=Vp;
            uz=Vz;
        elseif iaxe == 3
            up=Wp;
            uz=Wz;
        end
    else
        up=Up;
        uz=Uz;
    end
   %% Calcul des fonctions de Green
    makeglitch4 % la fonction de green est calculée dedans... 

%    save green_PL.mat green 
    ND2=15*floor(0.5/dt);
    NDsum=6*ND2;
    NDS=9*ND2;

    [am,im]=max(green); % recherche de la position max de la fonction de green

    isource=2;                     % s(t)= a*g(t) + b*(g'(t)) 

    if icorrect == 0
        if isource == 2
            NG1=2;            % instead of 3 or 4
            NG2=4;
            idelay=0;
        elseif isource == 3
            NG1=3;            % instead of 3 or 4
            NG2=5;
            idelay=1;
        end
    else
        NG1=4;                % instead of 3
        NG2=7;
    end

    if ifilter == 1
        for is=1:istep
            green=green-filtfilt(B,A,green);
            green_1=green_1-filtfilt(B,A,green_1);
            green_3=green_3-filtfilt(B,A,green_3);
            green_4=green_4-filtfilt(B,A,green_4);
            green_5=green_5-filtfilt(B,A,green_5);
        end
    elseif ifilter == 2
        for is=1:istep
            green=filtfilt(B,A,green);
            green_1=filtfilt(B,A,green_1);
            green_3=filtfilt(B,A,green_3);
            green_4=filtfilt(B,A,green_4);
            green_5=filtfilt(B,A,green_5);
        end
    end
    %%%%% for precursor

    % decoupage de la fonction de Green autour du max
    Green(:,1,iaxe)=green(im-ND2:im+NDS)/am;  % normalisation avec le maximum
                                              % On garde 15 points avant et 150
                                              % après.
    if NG1 == 2
        Green(:,2,iaxe)=green_1(im-ND2:im+NDS)/am;  % derivée de la fonction de Green
    %%%%%
    elseif NG1== 3 
        NGr=length(Green(:,1,iaxe));
        Green(1,2,iaxe)=Green(1,1,iaxe);
        Green(2:NGr,2,iaxe)=Green(1:NGr-1,1,iaxe);
        %%%%%
        Green(NGr,3,iaxe)=Green(NGr,1,iaxe);
        Green(1:NGr-1,3,iaxe)=Green(2:NGr,1,iaxe);
    end
    
    Weight(:,iaxe)=weight(im-ND2:im+NDS);

    ND=length(Green(:,1,iaxe));                        % longueur de le fonction de Green


    %weight(1:ND,iaxe)=1.;
    %weight(1:2*ND2+1,iaxe)=5.;

    % correct the green function 

    if icorrect == 1
        if dt == .5 
            load('dirstack/stacked_correction.mat');
        else
            load('dirstack20sps/stacked_correction.mat');
        end

    %%%% derivative is not efficient
        if length(stack_all0(:,iaxe)) == ND
            Green(:, NG1+1, iaxe) = stack_all0(:,iaxe);
        end
        NG1=NG1+1;
        NG2=NG2+1;

    end
    
    Green(1:ND, NG1+1, iaxe)=1.;
    Green(1:ND, NG1+2, iaxe)=[0:ND-1]*1;
    Green(1:ND, NG1+3, iaxe)=power([0:ND-1]*1, 2);
    
end % end of first loop on axes 


%%  start of inversion

%%%%% built the inverse problem
NN=length(u);
Na=(NG2+2);
i1=1; i2=ND;

Na3=Na*3;

%% initialize the inversion matrix

[aa, a]=initinversion(NG2,NDsum,Na,i1,i2,Green,Weight,recomp,epsi);
% aa est l'inverse de a

%%%%%%%%%%%% end of second loop on axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Var(1:NN,1:3)=1.1;
Var0(1:NN,1:3)=1.1;
Var1(1:NN,1:3)=.0;
gr=Green(1:NDsum, 1, 1);
rampg=gr(1)+(gr(NDsum)-gr(1))/(NDsum-1)*[0:NDsum-1];
im0=NDsum;

%%%%%% scann all selected max and invert for all selected max
for kk=NLOC:-1:1        % Boucle par amplitude décroissante des glitches
%%%% i is the starting point
	i=locm(kk)-ND2;
	if i >= 1 & i <= NN-ND
%%%%%%%%%%%%%%% loop for all glitches %%%%%%%%%%%%%%%%%%%%%%%
%%%%% check if i is the best
        rampu=u(i)+(u(i+NDsum-1)-u(i))/(NDsum-1)*[0:NDsum-1];
        rampv=v(i)+(v(i+NDsum-1)-v(i))/(NDsum-1)*[0:NDsum-1];
        rampw=w(i)+(w(i+NDsum-1)-w(i))/(NDsum-1)*[0:NDsum-1];
        Cu=xcorr(gr-rampg',u(i:i+NDsum-1)-rampu');
        Cv=xcorr(gr-rampg',v(i:i+NDsum-1)-rampv');
        Cw=xcorr(gr-rampg',w(i:i+NDsum-1)-rampw');
        [am1 imu]=max(abs(Cu));
        [am2 imv]=max(abs(Cv));
        [am3 imw]=max(abs(Cw));
        if abs(imu-imv) <= 2 & abs(imu-imw) > 2
                ima=floor((imu+imv)/2);
        elseif abs(imu-imw) <= 2 & abs(imu-imv) > 2
                ima=floor((imu+imw)/2);
        elseif abs(imv-imw) <= 2 & abs(imv-imu) > 2
                ima=floor((imv+imw)/2);
        else
                ima=floor((imu+imv+imw)/3);
        end
        Ishift(kk)=floor(ima-im0);
        if abs(Ishift(kk)) < 6 & i > Ishift(kk)
                i=i-Ishift(kk);
        end
%%%% recheck
	if i >= 1 & i <= NN-ND

		[cc,b]=inversion(NDsum, NG2, Na, i1, i2, Green, ...
                                        Weight, recomp, epsi, aa, u, v, w, uuz, i);
            % Y'a quoi dans b ??
            % CC => fonction  de green fittée sur le signal ? 
        
 		g0(1:NG2,i,1)=cc(1:NG2);
 		g0(1:NG2,i,2)=cc(Na+1:Na+NG2);
 		g0(1:NG2,i,3)=cc(2*Na+1:2*Na+NG2);

%%%%%%%%%%% END OF INVERSION %%%%%%%%%%%%

%%%%%%%%%%%% third loop on axes %%%%%%%%%%%%
		for iaxe=1:3

            %%%%%% compute synthetics glitches and drift
			synt(1:ND)=0.;
			Drift(1:ND)=0.;
            %%%%%%  synth is the synthetic with drift
			for j1=1:NG2
				synt(1:ND)=synt(1:ND)+cc(j1+(iaxe-1)*Na)*Green(1:ND,j1,iaxe)';
            end
            %%%%%  Drift is only the drift
			for j1=NG1+1:NG2
				Drift(1:ND)=Drift(1:ND)+cc(j1+(iaxe-1)*Na)*Green(1:ND,j1,iaxe)';
            end
            %%%% Delay through cross correlation
%		   [C21,lag21]=xcorr(synth(1:ND)-Drift(1:ND),s(i-1:i+ND-1,iaxe)-Drift(1:ND));

            % compute variance reduction
			var=0.;
			var0=0.;
			maxvar=0.;
			maxvar0=0.;

			for j=1:ND
				maxvar=max([maxvar,abs(s(i+j-1,iaxe)-synt(j))]);
				maxvar0=max([maxvar0,abs(s(i+j-1,iaxe)-Drift(j))]);
				var=var+Weight(j,iaxe)*(s(i+j-1,iaxe)-synt(j))^2;    % variance du signal sans drift et sans glitch
				var0=var0+Weight(j,iaxe)*(s(i+j-1,iaxe)-Drift(j))^2; % variance du signal sans le drift
            end
            
            
			Var(i,iaxe)=var/var0;  % réduction de variance pour chaque glitch.
			Var1(i,iaxe)=maxvar/maxvar0;
			Var0(i,iaxe)=var0;

            % check if results are model compatible


            %			if abs(g0(2,i,iaxe)/g0(1,i,iaxe)) > dtthreshold   % more than dtthreshold of offset: reject the glitch
            %				Var(i,iaxe)=1.1;
            %			end

			synt3D(ND2+1:ND+ND2,iaxe)=synt(1:ND);
			Drift3D(ND2+1:ND+ND2,iaxe)=Drift(1:ND);
			synt3D(1:ND2,iaxe)=NaN;
			Drift3D(1:ND2,iaxe)=NaN;

%%%%%%%%%%%% end of third loop on axes %%%%%%%%%%%%%%%%%%%%%%%%%
		end
%%%%%%%%%%% end loop for all glitches %%%%%%%%%%%%%%%%%%
	end
	end
end

%save deglitch7_g0.mat g0

Var0=Var;

%%%%% get the statistics and rank


%%%%% [unew,uglitch,udrift,udiff,u,uNGLITCH,uVar0,uVar1,val,val0,idx,Amp,uGreen,uND,ug0,ustack,ustack0,ND2,NDS,nstack,NG1,NG2]
%%%%%                                                              

%%%%%%%%%%%% forth loop on axes %%%%%%%%%%%

NNG=NN-ND;
nstack(1:3)=0;
%search=.5;
search = 1.;
%%%% rank the  glitch with the variance 

for iaxe=1:3

	ig=1;

	[val(1,iaxe),idd1] = min(Var(:,iaxe));   
	idx(1,iaxe)=idd1;
	val0(1,iaxe)=Var0(idd1,iaxe);
	ist=max([1,idd1-floor(ND/2)]);
	ist1=min([NN,idd1+floor(ND/2)]);

	Var(ist:ist1,iaxe) = search ;
    
    Amp(1:NG2, 1, iaxe)=g0(1:NG2, idx(1, iaxe), iaxe);  
    
    if idelay == 0
		Delay(ig,iaxe)=g0(2,idx(1,iaxe),iaxe)/g0(1,idx(ig,iaxe),iaxe);        %b/a
	elseif idelay == 1
		Delay(ig,iaxe)=dt*(g0(2,idx(ig,iaxe),iaxe)-g0(3,idx(1,iaxe),iaxe))/(g0(1,idx(ig,iaxe),iaxe)+g0(2,idx(ig,iaxe),iaxe)+g0(3,idx(ig,iaxe),iaxe));
	end

%%% go through all detected zero up to search threshold

	for i=2:length(locm)
  		[valval,idxidx] = min(Var(:,iaxe));
        %%% check if this is not too close from the previous glitchs
		if valval < search
			ig=ig+1;
			idx(ig,iaxe)=idxidx;
			val(ig,iaxe)=valval; % 
			val0(ig,iaxe)=Var0(idxidx,iaxe);
            % remove for the next iteration the time period of the glitch:
			ist=max([1,idxidx-floor(ND/2)]);
			ist1=min([NN,idxidx+floor(ND/2)]);
            %%%%% cas precedent
  			Var(ist:ist1,iaxe) = search;
  			Amp(1:NG2,ig,iaxe)=g0(1:NG2,idxidx,iaxe);
			if idelay == 0
				Delay(ig,iaxe)=g0(2,idxidx,iaxe)/g0(1,idxidx,iaxe);
			elseif idelay == 1
				Delay(ig,iaxe)=dt*(g0(2,idxidx,iaxe)-g0(3,idxidx,iaxe))/(g0(1,idxidx,iaxe)+g0(2,idxidx,iaxe)+g0(3,idxidx,iaxe));
			end
		end
    end %for i=2:length(locm) 

NZEROALL(iaxe)=ig; % ig: nombre de glitches

end

%NZERO

%save output_dgl7_B4_Ordering_by_var.mat Var0 locm Delay idx ND2 val Var1 Green Amp Var

%%%% alternative
Nmax=length(Var0(:,1));
Best(1:Nmax)=1.1;
Worst(1:Nmax)=1.1;
for i=1:Nmax
    Best(i)=min([Best(i),Var0(i,1),Var0(i,2),Var0(i,3)]);
    Worst(i)=max([Var0(i,1),Var0(i,2),Var0(i,3)]);
end

% save test2.mat Best
	
for iaxe=1:3

%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%

	snew(1:NN,iaxe)=0.;
	sglitch(1:NN,iaxe)=0.;
	sprecursor(1:NN,iaxe)=0;
	sdrift(1:NN,iaxe)=0.;
	sstack(1:ND,iaxe)=0.;
	sstack0(1:ND,iaxe)=0.;
	sstackp(1:ND,iaxe)=0.;
	sstack_res(1:ND,iaxe)=0.;
	jj=0;
	jjfirst=0;
	Time(1:ND)=[0:ND-1]*dt;

	inumber=0;
	gmax=max(abs(Green(1:ND,1,iaxe)));
	pmax=max(abs(green2sps_7(im-ND2:im+NDS)));
	corrg=3.71*cos(pi/180*30.);
	scale=10.;
	for i=1:NZEROALL(iaxe)
        LineDes=1;
        %%%%%% remove if only one of the the 3 axis is below threshold
    	%if Best(idx(i,iaxe)) < valthreshold  & val(i,iaxe) < 0.80
        %valthreshold = 0.05  ;                                                      % modif - Greg 6/10/2020
        if Best(idx(i,iaxe)) < valthreshold  & val(i,iaxe) <= 1.
%   		if Best(idx(i,iaxe)) < valthreshold  & val(i,iaxe) < 0.04
%   		if Best(idx(i,iaxe)) < valthreshold
            %%%%%% remove if this axis is below threshold
%   		if val(i,iaxe) < valthreshold 
%  		if Worst(i) < valthreshold 
			inumber=inumber+1;
			snew0(1:ND)=0.;
			sdrift0(1:ND)=0.;
			jj=jj+1;
			if jjfirst == 0
				LineDes=4;
				jjfirst=1;
			end
			jdx(jj,iaxe)=idx(i,iaxe);
            
            for j1=1:NG1		% paramètre a et b
                snew0(1:ND)=snew0(1:ND)+Amp(j1,i,iaxe)*Green(1:ND,j1,iaxe)'; %recomposition du signal de glitch ND longueur de la fn de Green
                                                                                                                 % i est le numéro du glitches

            end
            %%%%% Green(:,1,iaxe) is the response to H(t)
			[smax imax]=max(abs(snew0(1:ND))); smax=snew0(imax);
			SUVWmax(jj,iaxe)=smax;
        	for j1=NG1+1:NG2		
				sdrift0(1:ND)=sdrift0(1:ND)+Amp(j1,i,iaxe)*Green(1:ND,j1,iaxe)';
			end
			sglitch(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)=sglitch(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)+snew0(1:ND)'+sdrift0(1:ND)';
			sdrift(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)=sdrift0(1:ND)'; % idx = index de départ pour le début du glitch -> ce que je cherche pour replacer les glitches dans le temps
%%%%%
			ssnew0=snew0(1:ND)/smax;
			[snew1]=shift(ssnew0,Delay(i,iaxe),dt);
			[sdrift1]=shift(sdrift0,Delay(i,iaxe),dt);
			sstack0(1:ND,iaxe)=sstack0(1:ND,iaxe)+snew1(1:ND)';
			iplot = 0;
            if iplot == 1
				figure(100)
				subplot(3,6,1+(iaxe-1)*6)
 				plot(Time,snew1(1:ND))
				hold on
				xlabel('Second')
				ylim([-.75 1.25])
 				xlim([0 ND*dt])
				title('Modeled glitch')
				figure(102+iaxe)
				subplot(1,4,1)
  				plot(Time(1:ND),igg(iaxe)+jj+scale*snew1(1:ND)','k')
				hold on
				xlabel('Second')
				xlim([0 ND*dt])
			end
   			temp=(s(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)- ...
                            snew(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)-...
                            sdrift(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe))/smax;
                        
			[temp1]=shift(temp',Delay(i,iaxe),dt);
			if idx(i,iaxe) > ND
				NDshift=ND;
			else
				NDshift=idx(i,iaxe)-1;
			end
 			tempb=s(idx(i,iaxe)-NDshift:idx(i,iaxe)-1,iaxe)/smax;
			[temp1b]=shift(tempb',Delay(i,iaxe),dt);
			MaxAmp(i,iaxe)=max(temp);
			MinAmp(i,iaxe)=min(temp);
			sstack(1:ND,iaxe)=sstack(1:ND,iaxe)+(s(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)-snew(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)-sdrift(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe))/smax;
			nstack(iaxe)=nstack(iaxe)+1;
			snew(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)=snew(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)+snew0(1:ND)';
			sstack_res(1:ND,iaxe)=sstack_res(1:ND,iaxe)+temp1'-snew1';
			if iplot == 1
				figure(100)
				subplot(3,6,2+(iaxe-1)*6)
  				plot(Time,temp)
				hold on
				xlabel('Second')
				title('Observed glitch')
				ylim([-.75 1.25])
				xlim([0 ND*dt])
				figure(102+iaxe)
				subplot(1,4,1)
   				plot(Time,igg(iaxe)+jj+scale*temp1,'r','LineWidth',LineDes)
				xlim([0 ND*dt])
				subplot(1,4,2)
				plot(Time,igg(iaxe)+jj+2*scale*(temp1-snew1(1:ND)),'k')
				hold on
				xlim([0 ND*dt])
				subplot(1,4,3)
				plot(Time,igg(iaxe)+jj+2*scale*(temp1-snew1(1:ND)),'k')
				hold on
				xlim([0 ND*dt])
				figure(100)
			end
%%%%%% extract the spike at the begining
%%%%% remove the delta
%                      NDshift Nombre de point avant le début de séquece
%                      ND  longueur des fonction de green
%                      delay Poisition du max de la reponse à un heaviside
%                      g4delay Position du max de la reponse à un delta'	
			DeltaN=ND2-delay+g4delay;
		       % delta max
			NDGL=5;
 			diff(NDshift+1:ND+NDshift)=temp'-ssnew0;
 			diff(1:NDshift)=tempb';
			length(diff);
			Timeb(1:ND+NDshift)=([0:ND+NDshift-1]-floor((ND+NDshift)/2)-1)*dt;
			reducT=.5;
			if iaxe == 1
                reducu=1e9;
                for ic=-2:2
                    %%%% inversion entre idx(i,iaxe)-NDGL et idx(i,iaxe)+NDGL
                    [udgDt,udgDDt,guD0t,guD1t,var0t,vart,reducut]=betershift3(green2sps_7,green2sps_8,diff,NDshift+1+ic,dt,NDGL,DeltaN,g4delay);

                    reduc(ic+3)=reducut;
                    number(ic+3)=inumber;
                    if reducut < reducu
                        udgD=udgDt;
                        udgDD=udgDDt;
                        guD0=guD0t;
                        guD1=guD1t;
                        var0=var0t;
                        var=vart;
                        reducu=reducut;
                        ic1=NDshift+1+ic+DeltaN;
                    end
                end
            %%%%% precursor 
            %%%%% first data 1 is idx(i,iaxe)-NDshift 
            %%%%% data ic1 is therefore idx(i,iaxe)-NDshift  est la position 1
                ici=idx(i,iaxe)-NDshift + ic1 -1 ;
                %	if reducu < reducT
				sprecursor(ici-NDGL:ici+NDGL,1)=sprecursor(ici-NDGL:ici+NDGL,1)+udgD(ic1-NDGL:ic1+NDGL)'*smax;
                %	end
                %%%%% response to delta(t)
                [sdmax imax]=max(abs(udgD)); sdmax=udgD(imax);
                reducA(jj,1)=reducu;
			elseif iaxe == 2
                reducv=1e9;
                for ic=-2:2
                    [vdgDt,vdgDDt,gvD0t,gvD1t,var0t,vart,reducvt]=betershift3(green2sps_7,green2sps_8,diff,NDshift+1+ic,dt,NDGL,DeltaN,g4delay);
                    reduc(ic+3)=reducvt;
                    number(ic+3)=inumber;
                    if reducvt < reducv
                        vdgD=vdgDt;vdgDD=vdgDDt;
                        gvD0=gvD0t;gvD1=gvD1t;
                        var0=var0t;
                        var=vart;
                        reducv=reducvt;
                        ic1=NDshift+1+ic+DeltaN;
                    end
                end
                ici=idx(i,iaxe)-NDshift + ic1 -1 ;
    %			if reducv < reducT
                 sprecursor(ici-NDGL:ici+NDGL,2)=sprecursor(ici-NDGL:ici+NDGL,2)+vdgD(ic1-NDGL:ic1+NDGL)'*smax;
    %			end
                  [sdmax imax]=max(abs(vdgD)); 
                  sdmax=vdgD(imax);
                  reducA(jj,2)=reducv;
			elseif iaxe == 3
                reducw=1e9;
                for ic=-2:2
                [wdgDt,wdgDDt,gwD0t,gwD1t,var0t,vart,reducwt]=betershift3(green2sps_7,green2sps_8,diff,NDshift+1+ic,dt,NDGL,DeltaN,g4delay);
                reduc(ic+3)=reducwt;
                number(ic+3)=inumber;
                if reducwt < reducw
                    wdgD=wdgDt;wdgDD=wdgDDt;gwD0=gwD0t;gwD1=gwD1t;var0=var0t;var=vart;reducw=reducwt;
                    ic1=NDshift+1+ic+DeltaN;
                end
                end
                ici=idx(i,iaxe)-NDshift + ic1 -1 ;
    %			if reducwt < reducT
                    sprecursor(ici-NDGL:ici+NDGL,3)=sprecursor(ici-NDGL:ici+NDGL,3)+wdgD(ic1-NDGL:ic1+NDGL)'*smax;
    %			end
                            [sdmax imax]=max(abs(wdgD)); sdmax=wdgD(imax);
                reducA(jj,3)=reducw;
			end
%%%%% plot the delta for u,v,w
    		tempPrecurs=sprecursor(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)/smax;
			[tempPrecurs1]=shift(tempPrecurs',Delay(i,iaxe),dt);
			
            
            
            
            if iplot == 1
				figure(102+iaxe)
				subplot(1,4,2)
   				plot(Time(1:2*NDGL+1+DeltaN),igg(iaxe)+jj+scale*2*tempPrecurs1(1:2*NDGL+1+DeltaN),'r')
				subplot(1,4,4)
   				plot(Time,igg(iaxe)+jj+scale*2*(temp1-snew1(1:ND)-tempPrecurs1(1:ND)),'k')
				hold on
				xlim([0 ND*dt])
				figure(100)
			end
%%%%%%%%%%
			sstackp(1:ND,iaxe)=sstackp(1:ND,iaxe)+(s(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)-snew(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)-sdrift(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe)-sprecursor(idx(i,iaxe):idx(i,iaxe)+ND-1,iaxe))/smax;
			AmpSpike(jj,iaxe)=corrg*(sdmax/pmax)/(smax/gmax);
			subplot(3,6,5+(iaxe-1)*6)
			if sdmax*smax > 0 
                if reducA(jj,iaxe) < 0.25
                    DesCas='ro';
                elseif reducA(jj,iaxe) >= 0.25 & reducA(jj,iaxe) < 0.5
                    DesCas='r*';
                else
                    DesCas='r.';
                end
			else
                if reducA(jj,iaxe) < 0.25
                    DesCas='ko';
                elseif reducA(iaxe) >= 0.25 & reducA(jj,iaxe) < 0.5
                    DesCas='k*';
                else
                    DesCas='k.';
                end
			end
			if iplot == 1
                figure(100)
                semilogy(jj,abs(AmpSpike(jj,iaxe)),DesCas)
                hold on
                izap=1;
                if izap == 0
                figure(101)
                if iaxe == 1
                    subplot(4,1,1)
                    semilogy(24*(utlmst(idx(i,1))-iday),abs(AmpSpike(jj,1)),DesCas)
                    xlim([0 24])
                    xlabel('Local Time hr')
                    ylabel('Amplitude m')
                    title('Z or U component')
                    hold on
                elseif iaxe == 2
                    subplot(4,1,2)
                    semilogy(24*(vtlmst(idx(i,2))-iday),abs(AmpSpike(jj,2)),DesCas)
                    xlim([0 24])
                    xlabel('Local Time hr')
                    ylabel('Amplitude m')
                    title('X or V component')
                    hold on
                elseif iaxe == 3
                    subplot(4,1,3)
                    semilogy(24*(wtlmst(idx(i,3))-iday),abs(AmpSpike(jj,3)),DesCas)
                    xlim([0 24])
                    xlabel('Local Time hr')
                    ylabel('Amplitude m')
                    title('Y or W component')
                    hold on
                     %%%%
                    subplot(4,1,4)
                    if SUVWmax(jj,1) > 0 & jj > 0 & idx(i,1) > 0 & idx(i,1) <= length(utlmst)
                    semilogy(24*(utlmst(idx(i,1))-iday),SUVWmax(jj,1),'>r')
                    else jj > 0 & idx(i,1) > 0 & idx(i,1) <= length(utlmst)
                    semilogy(24*(utlmst(idx(i,1))-iday),-SUVWmax(jj,1),'<r')
                    end
                    hold on
                    if SUVWmax(jj,2) > 0  & jj > 0 & idx(i,2) > 0 & idx(i,2) <= length(vtlmst)
                    semilogy(24*(vtlmst(idx(i,2))-iday),SUVWmax(jj,2),'>g')
                    else jj > 0 & idx(i,2) > 0 & idx(i,2) <= length(vtlmst)
                    semilogy(24*(vtlmst(idx(i,2))-iday),-SUVWmax(jj,2),'<g')
                    end
                    if SUVWmax(jj,3) > 0 & jj > 0 & idx(i,3) > 0 & idx(i,3) <= length(wtlmst)
                    semilogy(24*(wtlmst(idx(i,3))-iday),SUVWmax(jj,3),'>b')
                    else jj > 0 & idx(i,3) > 0 & idx(i,3) <= length(wtlmst)
                    semilogy(24*(wtlmst(idx(i,3))-iday),-SUVWmax(jj,3),'<b')
                    end
                    xlim([0 24])
                end
                end
			end
			if iaxe == 3
%%%%%% comparing the orientations
				SUVWmax(jj,4)=sqrt(SUVWmax(jj,2)^2+SUVWmax(jj,3)^2);
				SUVWmax(jj,5)=sqrt(SUVWmax(jj,1)^2+SUVWmax(jj,2)^2+SUVWmax(jj,3)^2);
				AmpSpike(jj,4)=sqrt(AmpSpike(jj,2)^2+AmpSpike(jj,3)^2);
				AmpSpike(jj,5)=sqrt(AmpSpike(jj,1)^2+AmpSpike(jj,2)^2+AmpSpike(jj,3)^2);
				cos1(i)=(AmpSpike(jj,2)*SUVWmax(jj,2)+AmpSpike(jj,3)*SUVWmax(jj,3))/AmpSpike(jj,4)/SUVWmax(jj,4);
				cos2(i)=(AmpSpike(jj,1)*SUVWmax(jj,1)+AmpSpike(jj,2)*SUVWmax(jj,2)+AmpSpike(jj,3)*SUVWmax(jj,3))/AmpSpike(jj,5)/SUVWmax(jj,5);

			end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if iplot == 1
                figure(100)
                subplot(3,6,6+(iaxe-1)*6)
                plot(number,reduc,'k*')
                hold on
                plot(number(3),reduc(3),'rO')
                subplot(3,6,4+(iaxe-1)*6)
           
                NDp=length(diff);
                plot(Timeb(1:NDp)-DeltaN*dt+Delay(i,iaxe),diff(1:NDp),'k')
                hold on
                if iaxe == 1
                plot(Timeb(ic1-NDGL:ic1+NDGL)-DeltaN*dt+Delay(i,iaxe),udgDD(ic1-NDGL:ic1+NDGL),'r')
                elseif iaxe == 2
                plot(Timeb(ic1-NDGL:ic1+NDGL)-DeltaN*dt+Delay(i,iaxe),vdgDD(ic1-NDGL:ic1+NDGL),'r')
                elseif iaxe == 3
                plot(Timeb(ic1-NDGL:ic1+NDGL)-DeltaN*dt+Delay(i,iaxe),wdgDD(ic1-NDGL:ic1+NDGL),'r')
                end
                xlim([-5 5])
                ylim([-.2 .2])
                hold on
			end
			
        else
            %%%% reject the glitch
			Amp(:,i,iaxe)=0;   % modification par Greg 29/09/20 7:25
			val(i,iaxe)=0;
            %Amp(:,i,iaxe)=NaN(1);
			%val(i,iaxe)=NaN(1);
		end

	end
%%%%% number of glitvch
	NGLITCH(iaxe)=jj;

%%%%%% plot results
%	if iaxe == 3
%	if istep == 1
%	 	save totoAng1.mat SUVWmax AmpSpike cos1 cos2 reducA NGLITCH jdx
%	else
%		save totoAng2.mat SUVWmax AmpSpike cos1 cos2 reducA NGLITCH jdx
%	end
%	end
%%%%%%%


	if iplot == 1
		figure(100)
		subplot(3,6,3+(iaxe-1)*6)
		plot(sstack(:,iaxe)/nstack(iaxe),'k-','LineWidth',2)
		hold on
		plot(sstack0(:,iaxe)/nstack(iaxe),'r-','LineWidth',2)
		plot(5*sstackp(:,iaxe)/nstack(iaxe),'m-','LineWidth',2)
		plot(5*sstack_res(:,iaxe)/nstack(iaxe),'b-','LineWidth',1)
		xlabel('Second')
		title('Stacked glitch')
		legend('Data','Model','5xResidual')
		ylim([-.75 1.25])
	end
	nstack(iaxe)

%%%%
	sday=num2str(iday);
	sthres=num2str(valthreshold);
	sstep=num2str(istep);
	if icorrect == 0
		if dt == .5
			save([DIRSTACKF,'/stacked',sthres,'_',sstep,'sol',sday,'.mat'],'sstack','sstack0','nstack','sstack_res','valthreshold','dtthreshold','tfthreshold')
		elseif dt == 0.05
			save([DIRSTACK20SPS,'/stacked',sthres,'_',sstep,'sol',sday,'.mat'],'sstack','sstack0','nstack','sstack_res')
		end
	elseif icorrect == 1
		save([DIRSTACK_CORRECT_F,'/stacked',sthres,'_',sstep,'sol',sday,'.mat'],'sstack','sstack0','nstack','sstack_res','valthreshold','dtthreshold','tfthreshold')
	end

%save stack.mat ustack ustack0 nstack uz up ND2 NDS icase

	for i=1:length(u)
		if sglitch(i,iaxe) == 0
			sglitch(i,iaxe)=NaN(1);
			sdrift(i,iaxe)=NaN(1);
		end
	end

	if iaxe == 1
		unew=snew(:,1);
		uglitch=sglitch(:,1);
		udrift=sdrift(:,1);
		ustack=sstack(:,1);
		ustack0=sstack0(:,1);
		u=s(:,1);
	elseif iaxe == 2
		vnew=snew(:,2);
		vglitch=sglitch(:,2);
		vdrift=sdrift(:,2);
		vstack=sstack(:,2);
		vstack0=sstack0(:,2);
		v=s(:,2);
	elseif iaxe == 3
		wnew=snew(:,3);
		wglitch=sglitch(:,3);
		wdrift=sdrift(:,3);
		wstack=sstack(:,3);
		wstack0=sstack0(:,3);
		w=s(:,3);
%%%%%
		uuznew=recomp(1,1)*snew(:,1)+recomp(1,2)*snew(:,2)+recomp(1,3)*snew(:,3);
	end

end

%save recupData_deglitch7.mat snew s sglitch sdrift Best val sprecursor
%disp("File saved...")
end

