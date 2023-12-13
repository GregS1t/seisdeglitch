function [recomp] = rotate2zne_mat(dip1, azi1, dip2, azi2, dip3, azi3)


Geom= [sind(dip1) cosd(azi1)*cosd(dip1) cosd(dip1)*sind(azi1);
                sind(dip2) cosd(azi2)*cosd(dip2) cosd(dip2)*sind(azi2);
                sind(dip3) cosd(azi3)*cosd(dip3) cosd(dip3)*sind(azi3);
];
recomp=inv(Geom);
end