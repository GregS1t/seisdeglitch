function s2 = decimeFPGA(s1,ND,PFO,gain);
N1=length(s1);
N2=length(PFO);
ii=1;
for i=1:ND:N1
if i < N2
	s2(ii)=0.;
else
	s2(ii)=0.;
	for j=1:N2
		s2(ii)=s2(ii)+PFO(j)*s1(i-N2+j);
	end
end
s2(ii)=s2(ii)/gain;
ii=ii+1;
end
end

