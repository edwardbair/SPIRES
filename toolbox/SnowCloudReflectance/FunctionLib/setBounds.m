function [x0,lb,ub,variables] = setBounds(variables,fscript)
% set bounds for the unknowns

% all categories
snow = categorical({'snow'});
iceCloud = categorical({'iceCloud'});
waterCloud = categorical({'waterCloud'});
mixedCloud = categorical({'mixedCloud'});

S = SnowCloudLimits;

% if fSCA is a variable, and if more than one other endmember, then we need
% an additional unknown
if any(contains(variables,'fSCA','IgnoreCase',true))
    if size(fscript.R0,2)>1
        variables{end+1} = 'fEndMem1';
    end
end
x0 = zeros(1,length(variables));
lb = zeros(size(x0));
ub = zeros(size(x0));
for k=1:length(x0)
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
        case 'dust'
            x0(k) = mean(S.dust);
            lb(k) = S.dust(1);
            ub(k) = S.dust(2);
        case 'dustRadius'
            x0(k) = S.defaultDustRadius;
            lb(k) = S.dustRadius(1);
            ub(k) = S.dustRadius(2);
        case 'soot'
            x0(k) = mean(S.soot);
            lb(k) = S.soot(1);
            ub(k) = S.soot(2);
        case 'sootRadius'
            x0(k) = S.defaultSootRadius;
            lb(k) = S.sootRadius(1);
            ub(k) = S.sootRadius(2);
        case 'fSCA'
            x0(k) = 0.5;
            lb(k) = 0;
            ub(k) = 1;
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
        case 'nRwt'
            lb(k) = 0;
            ub(k) = 1;
            x0(k) = 0.5;
        otherwise
            if contains(variables{k},'fEndMem','IgnoreCase',true)
                x0(k) = 0.5;
                lb(k) = 0;
                ub(k) = 1;
            else
                error('variable ''%s'' not recognized for setting bounds',variables{k})
            end
    end
end