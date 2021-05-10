load('BeresData.mat')
SZ = 26.74;
t = BeresData.cleanRefl>0 & BeresData.waveL_um>=0.9; % get rid of negative refl and use wavelengths only beyone 0.9 um
% clean snow, but probably contains dust
disp('"clean" snow, probably contains dust but solving for grain size only')
[oSclean,statsClean,Pclean] = invertSnowCloudSpectralRefl(BeresData.cleanRefl(t),{'ssa'},'snow','wavel',BeresData.waveL_um(t),'waveu','um','cosZ',cosd(SZ) );
disp(oSclean)
disp(statsClean)
effective_radius_clean = SSA2radius(oSclean.ssa,'um') %#ok<NOPTS>
% sooty snow
disp('snow with BC, along with dust')
tbc = BeresData.BCRefl>0; % get rid of negative refl, but include whole spectrum
[oSBC, statsBC, PBC] = invertSnowCloudSpectralRefl(BeresData.BCRefl(tbc),{'soot','dust'},'snow','ssa',oSclean.ssa,'wavel',BeresData.waveL_um(tbc),'waveu','um','cosZ',cosd(SZ) );
disp(oSBC)
disp(statsBC)
disp('snow with BC, assuming no dust')
[oSBCnd,statsBConly] = invertSnowCloudSpectralRefl(BeresData.BCRefl(tbc),{'soot'},'snow','ssa',oSclean.ssa,'wavel',BeresData.waveL_um(tbc),'waveu','um','cosZ',cosd(SZ) );
disp(oSBCnd)
disp(statsBConly)

% plots
Rclean = SnowCloudSpectralRefl(Pclean,'wavelength',BeresData.waveL_um,'waveu','um');
plot(BeresData.waveL_um(t),BeresData.cleanRefl(t),BeresData.waveL_um,Rclean.refl)
Rsoot = SnowCloudSpectralRefl(PBC,'wavelength',BeresData.waveL_um,'waveu','um');
hold on;
plot(BeresData.waveL_um(tbc),BeresData.BCRefl(tbc),BeresData.waveL_um,Rsoot.refl)
legend({'measClean','modelClean','measBC','modelBC'})

