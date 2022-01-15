function [newU,newP] = fixShade(unknowns,Prescription)
%get rid of a shade endmember by adding a new R0 and adding a new fraction
%to the fSCA vector

newU = unknowns;
newP = Prescription;
if Prescription.substance == categorical({'snow'})
    t = contains(unknowns,'shade','IgnoreCase',true);
    if any(t)
        newU = unknowns(~t);
        x = [min(newP.Spectrum.wavelength) max(newP.Spectrum.wavelength)];
        y = [0 0];
        F = fit(x',y','poly1');
        newP.Substrate.fR0 = cat(2,newP.Substrate.fR0,F);
        newP.snow.fSCA = [newP.snow.fSCA(1) 1-sum(newP.snow.fSCA)];
    end
end
end