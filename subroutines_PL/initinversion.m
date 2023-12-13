function [aa, a, i1, i2]=initinversion(NG2,NDsum,Na,  i1, i2, Green, Weight, recomp, epsi)

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
% --------
% OUTPUT
% 


Na3=3*Na;

a(1:Na3,1:Na3)=0.;

%%%%%%%%%%% cross terms for vertical axis -> 
for iaxe=1:3
	for jaxe=1:3
        for j2=1:NG2
            temp(1:NDsum)=recomp(1, jaxe)*Weight(1:NDsum, jaxe).*Green(1:NDsum, j2, jaxe);
            for j1=1:NG2
                a(j1+(iaxe-1)*Na, j2+(jaxe-1)*Na)=epsi*recomp(1,iaxe)*dot(Green(1:NDsum,j1,iaxe),temp(1:NDsum));
            end
        end
	end
end


%%%%%%%%%%%% second loop on axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Lagrange multiplier to lock start and end points

for iaxe=1:3
	for j2=1:NG2
		temp(1:NDsum)=Weight(1:NDsum, iaxe).*Green(1:NDsum, j2, iaxe);
		for j1=1:NG2
			a(j1+(iaxe-1)*Na,j2+(iaxe-1)*Na)=a(j1+(iaxe-1)*Na,j2+(iaxe-1)*Na)+dot(Green(1:NDsum,j1,iaxe), temp(1:NDsum));
		end
	end
%%%%%
	for j1=1:NG2
		a(NG2+1+(iaxe-1)*Na, j1+(iaxe-1)*Na)=Green(i1, j1,iaxe)/2; 
        a(j1+(iaxe-1)*Na, NG2+1+(iaxe-1)*Na)=a(NG2+1+(iaxe-1)*Na, j1+(iaxe-1)*Na);              % impose le debut de glitch identique
		
        a(NG2+2+(iaxe-1)*Na,j1+(iaxe-1)*Na)=Green(i2, j1,iaxe)/2; 
        a(j1+(iaxe-1)*Na, NG2+2+(iaxe-1)*Na)=a(NG2+2+(iaxe-1)*Na, j1+(iaxe-1)*Na);
	end
end
%%%%%%%%%%%% end of second loop on axes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aa=inv(a(1:Na3, 1:Na3));
%save a_initinversion.mat a aa
end
