function [a] = shift(b,delta,dt)
	spec=fft(b);
	N=length(b);
	NF=floor(N/2)+1;
	df=1/(N*dt);
	freq=2*pi*[0:NF-1]*df*1i*delta;
	spec(1:NF)=spec(1:NF).*exp(-freq);
	a=ifft(spec,'symmetric');
	end
