function [ts_c2] = remove_tfvit(ts_or,t,dt,gain,z,p,A0,fc)
%this removes the instrument response

%%%Find gaps
N1=[];
N2=[];
N1=1;
tol=.5;
for ii=2:length(t)
    if 86400*(t(ii)-t(ii-1))>((1+tol)*dt)
        N2=[N2;ii-1];
        N1=[N1;ii];
    end
end
    N2=[N2; length(t)];
N1;
N2;
NGAPS=length(N2); %Number of gaps+1


%%%Loop on data without gaps
ts_c2=[];
for NT=1:NGAPS
clear w tf ts
ts0=ts_or(N1(NT):N2(NT));
%%%%%% make the connection between first and last point
NS=length(ts0);
trend=ts0(1)+[1:NS]'*(ts0(NS)-ts0(1))/(NS-1);
ts(1:floor(NS/2))=trend(floor(NS/2):-1:1);
ts(floor(NS/2)+1:floor(NS/2)+NS)=ts0;
ts(floor(NS/2)+NS+1:2*NS)=trend(floor(NS/2)+1:NS);
NS1=NS;
NS=2*NS;
fft_ts=fft(ts);
fft_ts(1)=0;
%%%%%% modif par PL
NF=floor(NS/2+1);
for i=1:NF
w(i)=(i-1)/NS/dt*2*pi*1i;
tf(i)=1.;
end
tf0=1.;
w0=.1*2*pi*1i;

% p=p*(2*pi);
% z=z*(2*pi);
for ii=1:length(z)
    tf=tf.*(w-z(ii));
    tf0=tf0*(w0-z(ii));
end

for ii=1:length(p)
    tf=tf./(w-p(ii));
    tf0=tf0/(w0-p(ii));
end

%save toto.mat gain w tf A0 

tf=gain*A0*tf;

figure(55)

loglog(abs(w)/2/pi,abs(tf))
hold on
xlabel('Frequency')
ylabel('Transfer function DU/(m/s)')

%%%%%

for i=1:NF
fft_ts(i)=fft_ts(i)*w(i)/tf(i);
end

fft_ts(1)=0;

ts_c=ifft(fft_ts,'symmetric');

[B,A]=butter(2,fc*dt,'high');         %%%% high pass freq
ts_c=(filtfilt(B,A,ts_c));            %%%%

ts_c2=[ts_c2;ts_c(floor(NS1/2)+1:floor(NS1/2)+NS1)'];

end
end

