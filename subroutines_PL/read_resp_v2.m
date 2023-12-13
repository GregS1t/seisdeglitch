function [ gain z p units A0] = read_resp( filename )
%This script reads the response file and outputs gain, poles, zeros, and
%units 
%NO TIME or FIR filters are considered (only the Gain of the FIR)
fid=fopen(filename,'r');

a=1;
while(a)
   c=fgetl(fid);

   k=strfind(c,'Response in units');
   if length(k)>0
       a=0;
       kk=strfind(c,'DISPLACEMENT');
       if length(kk)>0
           units='DIS';
       end
       kk=strfind(c,'VELOCITY');
       if length(kk)>0
           units='VEL';
       end
       kk=strfind(c,'ACCELERATION');
       if length(kk)>0
           units='ACC';
       end
       kk=strfind(c,'COUNTS');
       if length(kk)>0
           units='CTS';
       end
   end
end
c=fgetl(fid);
c=fgetl(fid);
A0=str2num(c(52:end));
c=fgetl(fid);
c=fgetl(fid);
nz=str2num(c(52:end));
c=fgetl(fid);
np=str2num(c(52:end));

%%read zeros and poles
c=fgetl(fid);
c=fgetl(fid);
z=[];
for ii=1:nz
    c=fgetl(fid);
    zr=str2num(c(17:29));
    zi=str2num(c(31:43));
    z=[z zr+1i*zi];
end
c=fgetl(fid);
c=fgetl(fid);
p=[];
for ii=1:np
    c=fgetl(fid);
    pr=str2num(c(17:29));
    pi=str2num(c(31:43));
    p=[p pr+1i*pi];
end
a=1;
while(a)
    c=fgetl(fid);
    k=strfind(c,'Gain:');
    if length(k)>0
        a=0;
        gain_stage1=str2num(c(52:end));
    end
end
a=1;
while(a)
    c=fgetl(fid);
    k=strfind(c,'Gain:');
    if length(k)>0
        a=0;
        gain_stage2=str2num(c(52:end));
    end
end

gain=gain_stage1*gain_stage2;

a=1;

    while(a)
        if strcmp(class(c),'double')
            a=0;
        end
        c=fgetl(fid);
        k=strfind(c,'Start date');
        if length(k)>0
            a=0;
        end
        k=strfind(c,'Gain:');
        if length(k)>0
                fir_gain=str2num(c(52:end));
                gain=gain*fir_gain;
        end
    end
    
fclose(fid);


end

