function [ gain z p units A0] = read_resp( filename )
%This script reads the response file and outputs gain, poles, zeros, and
%units 
fid=fopen(filename,'r');
a=1;
while(a)
   c=fgetl(fid);
   
   k=strfind(c,'Response in units');
   if length(k)>0
%       a=0;
       kk=strfind(c,'DISPLACEMENT');
       if length(kk)>0
           units='DIS';
        a=0;
       end
       kk=strfind(c,'VELOCITY');
       if length(kk)>0
           units='VEL';
        a=0;
       end
       kk=strfind(c,'ACCELERATION');
       if length(kk)>0
           units='ACC';
        a=0;
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
    
end

