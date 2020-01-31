function postprocess_Horizons(demfile, whichBlock )
% postprocess_Horizons(demfile, whichBlock [,folder] )
%   Use the results from prep_Horizons to smooth and interpolate horizon
%   vector and calculate view factors
%
% Input
%   demfile - name of original demfile
%   whichBlock - which Hblock file, runs in parallel
%
% Output
%   written to files in same folder as demfile
%   each contains original horizon header and piece of horizon grid
%

doSmooth = true; % flag to turn on/off smoothing
[folder,demname,ext] = fileparts(demfile); %#ok<ASGLU>
if isempty(folder)
    folder = '.';
end

if nargin<2
    if isdeployed
        disp(['usage: ' mfilename ' demfile whichBlock'])
    else
        disp(['usage: ' mfilename '(demfile, whichBlock)'])
    end
    return
end

% which block to read
if ~isnumeric(whichBlock)
    whichBlock = str2double(whichBlock);
end

% output diary file
if isdeployed
    diaryFolder = fullfile(folder,'output');
    if exist(diaryFolder,'file')~=7
        mkdir(diaryFolder);
    end
    StartAzureDiary(mfilename,diaryFolder,demname,...
        ['_Hblock' num2str(whichBlock,'%04d')]);
end

f = dir(fullfile(folder,[demname '*' ['Hblock' num2str(whichBlock,'%04d')] '*.mat']));

assert(length(f)==1,'should be just one file with name %s',...
    [demname '*' ['Hblock' num2str(whichBlock,'%04d')] '*.mat'])
outputfile = [demname '_Hfilter' num2str(whichBlock,'%04d') '.mat'];
disp(['filtering file ' f.name])

% input file, and make sure we got the variables
H = load(fullfile(folder,f.name));
assert(isfield(H,'azm') && isfield(H,'sinH') &&...
    isfield(H,'aspect') && isfield(H,'slope'),...
    'check file %s, either azm, sinH, aspect, or slope does not exist',...
    f.name)

% conversion from floating point to integer, using int16 even though most
% values are unsigned to allow for flexibility
divisor = double(intmax('int16'));

% output horizon azimuths
xp = (-180:180)';
output = zeros(length(xp),size(H.sinH,2),'int16');
V = zeros(1,size(H.sinH,2));
tempS = H.slope;
tempA = H.aspect;

parfor k=1:size(H.sinH,2)
    pair = sortrows([H.azm(:,k) H.sinH(:,k)]); %#ok<PFBNS>
    x = pair(:,1);
    y = pair(:,2);
    [a,ia] = unique(x);
    if length(a)~=length(x)
        x = a;
        y = y(ia);
    end
    
    % intperpolate for every degree
    F = griddedInterpolant(x,y,'pchip','nearest');
    yp = F(xp);
    % smooth?
    if doSmooth
        yp = smooth(yp,5,'lowess');
    end
    % force periodic
    endpt = mean([yp(1) yp(end)]);
    yp(1) = endpt;
    yp(end) = endpt;
    
    yp(yp<0) = 0; % correct for undershoot
    V(k) = viewf(yp,tempS(k),tempA(k));
    output(:,k) = int16(round(yp*divisor));
end
Horz.sinH = output;
Horz.V = V;
B = whos('Horz');
if B.bytes>2^31
    save(fullfile(folder,outputfile),'-struct','Horz','-v7.3')
else
    save(fullfile(folder,outputfile),'-struct','Horz')
end
disp(['file ' outputfile ' saved'])

end