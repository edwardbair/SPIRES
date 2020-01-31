function Tbl = summarizeSnowCloud(cloud,snow)
%summarize results of comparisons of retrievals and cloud and snow
%properties from the same spectrum (i.e., from similarSnowCloud.m)
%
%Input
%   cloud - cloud structure from similarSnowCloud.m
%   snow - snow structure from similarSnowCloud.m
%

Tbl = [];
for c=1:size(cloud,2)
    for r=1:size(cloud,1)
        if isempty(cloud(r,c).outputStruct) || isempty(snow(r,c).outputStruct)
            continue;
        end
        tc = cloud(r,c);
        dt = datetime(tc.dnval,'ConvertFrom','datenum');
        dt.Format = 'dd-MMM-uuuu';
        thisRow = table({tc.file},dt,...
            [int16(str2double(tc.file(4:6))) int16(str2double(tc.file(7:9)))],...
            uint32(round(tc.idx)),{tc.ErrorType},tc.solarZ,...
            'VariableNames',...
            {'file','date','pathrow','pixelIDX','ErrorType','solarZ'});
        thisTbl = [];
        for s=1:length(cloud(r,c).outputStruct)
            oc = cloud(r,c).outputStruct(s);
            os = snow(r,c).outputStruct(s);
            if isempty(oc) || isempty(os)
                continue;
            end
            if length(os.oStruct.otherEndMem)==1
                os.oStruct.otherEndMem(2) = 0;
            end
            if isfield(oc.oStruct,'wetness')
                cloudResult = table({oc.Background},oc.stats.resnorm,...
                    oc.stats.exitflag,oc.oStruct.radius,...
                    oc.oStruct.waterEquivalent,oc.oStruct.wetness,...
                    'VariableNames',{'Background',...
                    'Cresnorm','Cexitflag','Cradius','waterEq','wetness'});
            else
                cloudResult = table({oc.Background},oc.stats.resnorm,...
                    oc.stats.exitflag,oc.oStruct.radius,...
                    oc.oStruct.waterEquivalent,...
                    'VariableNames',{'Background',...
                    'Cresnorm','Cexitflag','Cradius','waterEq'});
            end
            snowResult = table(os.stats.resnorm,...
                os.stats.exitflag,os.oStruct.radius,...
                os.oStruct.fSCA,os.oStruct.dust,os.oStruct.otherEndMem,...
                'VariableNames',{'Sresnorm','Sexitflag','Sradius','fSCA','dust','otherEndMem'});
            residTbl = table(cloud(r,c).refl,oc.stats.residual',os.stats.residual',...
                'VariableNames',{'measuredReflectance','cloudResidual','snowResidual'});
            % leave out those where snow radius and cloud radius are max'd
            if ~(oc.oStruct.radius>49 && os.oStruct.radius>1990)
                thisTbl = [thisTbl;...
                    [thisRow cloudResult snowResult residTbl]]; %#ok<AGROW>
            end
        end
        Tbl = [Tbl;thisTbl]; %#ok<AGROW>
        
    end
end
end