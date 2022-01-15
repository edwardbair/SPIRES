function [x0,lb,ub,variables] = setBounds(variables,fscript,solver)
% set bounds for the unknowns

% all categories
snow = categorical({'snow'});
iceCloud = categorical({'iceCloud'});
waterCloud = categorical({'waterCloud'});
mixedCloud = categorical({'mixedCloud'});

S = SnowCloudLimits;

% if fSCA is a variable, it's at the end (from setUnknowns), so increment x0
% and solve in a constrained way so that sum(all fractions)=1, but this
% requires using fmincon not lsqnonlin
% 0.5 is default for some variables
if any(contains(variables,'fSCA'))
    if strcmpi(solver,'leastSquares')
        x0 = 0.5+zeros(1,length(variables));
        assert(length(fscript.snow.fSCA)==2,...
            'can''t solve for multiple endmembers with ''leastSquares'' option, use ''normResiduals'' or ''spectralAngle'' instead')
    else
        x0 = 1/length(fscript.snow.fSCA)+zeros(1,length(variables)+length(fscript.snow.fSCA)-1);
    end
else
    x0 = zeros(1,length(variables));
end

% some variables have [0 1] lower and upper defaults, so set those
% initially
lb = zeros(size(x0));
ub = ones(size(x0));
for k=1:length(variables)
    switch variables{k}
        case 'radius'
            switch fscript.substance
                case snow
                    x0(k) = S.defaultSnowRadius;
                    lb(k) = S.snowRadius(1);
                    ub(k) = S.snowRadius(2);
                case {iceCloud,mixedCloud}
                    x0(k) = S.defaultIceCloudRadius;
                    lb(k) = S.iceCloudRadius(1);
                    ub(k) = S.iceCloudRadius(2);
                case waterCloud
                    x0(k) = S.defaultWaterCloudRadius;
                    lb(k) = S.waterCloudRadius(1);
                    ub(k) = S.waterCloudRadius(2);
            end
        case 'waterRadius'
            switch fscript.substance
                case mixedCloud
                    x0(k) = S.defaultWaterCloudRadius;
                    lb(k) = S.waterCloudRadius(1);
                    ub(k) = S.waterCloudRadius(2);
                case snow
                    x0(k) = S.defaultWaterInSnowRadius;
                    lb(k) = S.waterInSnowRadius(1);
                    ub(k) = S.waterInSnowRadius(2);
                otherwise
                    error('''waterRadius'' not an appropriate unknown for substance %s',...
                        fscript.substance)
            end
        case 'wetness'
            switch fscript.substance
                case snow
                    x0(k) = mean(S.wetSnow);
                    lb(k) = S.wetSnow(1);
                    ub(k) = S.wetSnow(2);
                otherwise
                    x0(k) = 0.5;
                    lb(k) = 0;
                    ub(k) = 1;
            end
        case 'LAPfraction'
            if strcmpi(fscript.snow.LAP,'dust')
                x0(k) = mean(S.dust);
                lb(k) = S.dust(1);
                ub(k) = S.dust(2);
            elseif strcmpi(fscript.snow.LAP,'soot')
                x0(k) = mean(S.soot);
                lb(k) = S.soot(1);
                ub(k) = S.soot(2);
            else % use range of dust, soot
                x0(k) = mean(S.dust);
                lb(k) = S.soot(1);
                ub(k) = S.dust(2);
            end
        case {'LAPfraction1','LAPfraction2'}
            if contains(variables{k},'1')
                nv = 1;
            else
                nv = 2;
            end
            if strcmpi(fscript.snow.LAP{nv},'dust')
                x0(k) = mean(S.dust);
                lb(k) = S.dust(1);
                ub(k) = S.dust(2);
            elseif strcmpi(fscript.snow.LAP{nv},'soot')
                x0(k) = mean(S.soot);
                lb(k) = S.soot(1);
                ub(k) = S.soot(2);
            else % use range of dust, soot
                x0(k) = mean(S.dust);
                lb(k) = S.soot(1);
                ub(k) = S.dust(2);
            end
        case 'LAPradius'
            if strcmpi(fscript.snow.LAP,'dust')
                x0(k) = S.defaultDustRadius;
                lb(k) = S.dustRadius(1);
                ub(k) = S.dustRadius(2);
            elseif strcmpi(fscript.snow.LAP,'soot')
                x0(k) = S.defaultSootRadius;
                lb(k) = S.sootRadius(1);
                ub(k) = S.sootRadius(2);
            else % use range of dust, soot
                x0(k) = S.defaultDustRadius;
                lb(k) = S.sootRadius(1);
                ub(k) = S.dustRadius(2);
            end
        case {'LAPradius1','LAPradius2'}
            if contains(variables{k},'1')
                nv = 1;
            else
                nv = 2;
            end
            if strcmpi(fscript.snow.LAP{nv},'dust')
                x0(k) = S.defaultDustRadius;
                lb(k) = S.dustRadius(1);
                ub(k) = S.dustRadius(2);
            elseif strcmpi(fscript.snow.LAP{nv},'soot')
                x0(k) = S.defaultSootRadius;
                lb(k) = S.sootRadius(1);
                ub(k) = S.sootRadius(2);
            else % use range of dust, soot
                x0(k) = S.defaultDustRadius;
                lb(k) = S.sootRadius(1);
                ub(k) = S.dustRadius(2);
            end
        case 'waterEquivalent'
            switch fscript.substance
                case snow
                    x0(k) = 100;
                    lb(k) = 5;
                    ub(k) = Inf;
                case {waterCloud,mixedCloud}
                    x0(k) = mean(S.waterCloudWE);
                    lb(k) = S.waterCloudWE(1);
                    ub(k) = S.waterCloudWE(2);
                case iceCloud
                    x0(k) = mean(S.iceCloudWE);
                    lb(k) = S.iceCloudWE(1);
                    ub(k) = S.iceCloudWE(2);
            end
        case 'muS'
            x0(k) = 0.5;
            % defaults for lb and ub
        case 'fSCA'
            % use defaults so don't need to set
        otherwise
            error('variable ''%s'' not recognized for setting bounds',variables{k})
    end
end
end