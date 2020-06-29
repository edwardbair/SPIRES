function vgf=GOvgf(cc,theta_s,phi_s,theta_v,phi_v,b_R)
%Liu et al 2004, Eq 8, DOI: 10.1002/hyp.5802
%calculate viewable gap fraction based on viewing geometry
%input, scalars or vectors or matrices: 
%cc: canopy cover fraction
%theta_s - slope angle, deg
%phi_s - slope aspect, deg
%theta_v - view zenith, deg
%phi_v -view azimuth, deg
% b_R - b/R - avg vertical crown radius divided by 
%average horizontal crown radius 

theta_v_prime=atand(b_R.*tand(theta_v));
theta_s_prime=90-atand((b_R.*tand(90-theta_s)));

phi_v_prime=phi_v-phi_s;
% theta_v_prime_prime=acosd(cosd(phi_v_prime.*sind(theta_v_prime).*...
%     sind(theta_s_prime)+cosd(theta_v_prime).*cosd(theta_s_prime));

vgf=(1-cc).^((cosd(theta_s_prime))./...
    (cosd(phi_v_prime).*sind(theta_v_prime).*sind(theta_s_prime)+...
    cosd(theta_v_prime).*cosd(theta_s_prime)));