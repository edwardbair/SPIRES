function scaledValues=bandScaling(sensor,bandNo,variable,inputValues)

% scale the values to make the interpolation quasi-linear but must keep the
% sorting progression, hence the use of the minus signs
if isnumeric(bandNo)
    bandNo = categorical(bandNo);
elseif ischar(bandNo)
    bandNo = categorical({bandNo});
end
switch lower(sensor)
    case 'landsatoli'
        switch lower(variable)
            case 'ssa'
                switch bandNo
                    case {categorical(1),categorical(2),categorical(3),...
                            categorical(4),categorical(5),categorical(8)}
                        scaledValues = -1./sqrt(inputValues);
                    case categorical(9)
                        scaledValues = inputValues.^0.01;
                    otherwise %linear
                        scaledValues = inputValues;
                end
            case {'dust','soot'}
                switch bandNo
                    case {categorical(1),categorical(2),categorical(3),...
                            categorical(8)}
                        scaledValues = sqrt(inputValues);
                    case categorical(4)
                        scaledValues = inputValues.^0.75;
                    otherwise %linear
                        scaledValues = inputValues;
                end
            otherwise
                % no change
                scaledValues = inputValues;
        end
    case 'modis'
        switch lower(variable)
            case 'ssa'
                switch bandNo
                    case {categorical(1),categorical(2),categorical(3),...
                            categorical(4)}
                        scaledValues = -1./sqrt(inputValues);
                    case categorical(5)
                        scaledValues = inputValues.^0.04;
                    otherwise %linear
                        scaledValues = inputValues;
                end
            case {'dust','soot'}
                switch bandNo
                    case {categorical(3),categorical(4)}
                        scaledValues = sqrt(inputValues);
                    case categorical(1)
                        scaledValues = inputValues.^0.7;
                    otherwise %linear
                        scaledValues = inputValues;
                end
            otherwise
                % no change
                scaledValues = inputValues;
        end
end
end