function [cc, b]=inversion(NDsum, NG2, Na, i1, i2, Green, Weight, recomp, epsi, aa, u, v, w, uuz, i)



%%%%%%%%%%%% start of inversion %%%%%%%%
% NG2      = 4 ou 5 selon le type de modélisation dans notre cas, = 2
% NDsum = 6*ND2;
% ND2      = 15*floor(0.5/dt);
% Na         = (NG2+2);               % s(t)= a*g(t) + b*(g'(t)) -> a, b, moyenne, pente -> 4 coefs 
% i1           = 1;
% i2           = ND (ND=length(Green(:,1,iaxe))
% Green    = Fonction de green calculée dans makeglitch4
% Weight   = Poids calculés autour de im (position du maximum de la fonction de Green)
% Recomp = matrice de rotation
% epsi        = 0 dans notre cas
%



Na3=3*Na;

		for j1=1:NG2
			temp(1:NDsum)=Weight(1:NDsum,1).*u(i:i+NDsum-1);
			tempZ(1:NDsum)=Weight(1:NDsum,1).*uuz(i:i+NDsum-1);
			b(j1)=dot(Green(1:NDsum,j1,1),temp(1:NDsum));
			b(j1)=b(j1)+epsi*recomp(1,1)*dot(Green(1:NDsum,j1,1),tempZ(1:NDsum));
%%%%%
			temp(1:NDsum)=Weight(1:NDsum,2).*v(i:i+NDsum-1);
			tempZ(1:NDsum)=Weight(1:NDsum,2).*uuz(i:i+NDsum-1);
			b(j1+Na)=dot(Green(1:NDsum,j1,2),temp(1:NDsum));
			b(j1+Na)=b(j1+Na)+epsi*recomp(1,2)*dot(Green(1:NDsum,j1,2),tempZ(1:NDsum));
%%%%%
			temp(1:NDsum)=Weight(1:NDsum,3).*w(i:i+NDsum-1);
			tempZ(1:NDsum)=Weight(1:NDsum,3).*uuz(i:i+NDsum-1);
			b(j1+2*Na)=dot(Green(1:NDsum,j1,3),temp(1:NDsum));
			b(j1+2*Na)=b(j1+2*Na)+epsi*recomp(1,3)*dot(Green(1:NDsum,j1,3),tempZ(1:NDsum));
		end
%%%%%% impose the start and end on the trend
 		b(NG2+1)=u(i+i1-1)/2;
 		b(NG2+2)=u(i+i2-1)/2;
 		b(NG2+1+Na)=v(i+i1-1)/2;
 		b(NG2+2+Na)=v(i+i2-1)/2;
 		b(NG2+1+2*Na)=w(i+i1-1)/2;
 		b(NG2+2+2*Na)=w(i+i2-1)/2;
		cc=(aa(1:Na3,1:Na3)*b(1:Na3)')';

end

