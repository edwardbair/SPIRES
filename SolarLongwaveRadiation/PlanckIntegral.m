function y=PlanckIntegral(a,b,T)
% y=PlanckIntegral(a, b, temperature)
% compute integral of Planck function from wavelengths a to b at temp T
% a and b are in meters
% temperature in Kelvin
% result is in W/m^3/sr

tval = T;

y=integral(@pfun,a,b,'ArrayValued',true);

    function s=pfun(w)
        s=planck(w,tval);
    end
end