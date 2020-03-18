function R = edgint( ABVEDG, BLWEDG, AEUP, BEUP,ACCUR)
% R = edgint( ABVEDG, BLWEDG, AEUP, BEUP)
%        Computes edge correction integrals to geometric optics
% Input
%	ABVEDG = Integrand function for above-edge integral
%	BLWEDG = Integrand function for below-edge integral
%	AEUP   = Upper limit for above-edge integral
%	BEUP   = Upper limit for below-edge integral

above = quadgk(ABVEDG,0,AEUP,'RelTol',ACCUR);
below = quadgk(BLWEDG,0,BEUP,'RelTol',ACCUR);

R = above+below;

end