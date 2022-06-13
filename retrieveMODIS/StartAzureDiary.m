function [ output_host ] = StartAzureDiary(inputFunction,outputLocation,varargin)
%create diary file for Azure instance, but only if running deployed
% [ output_host ] = StartAzureDiary(inputFunction,outputLocation,varargin)
%
% Input, 2 required, in this order
%   inputFunction - perhaps from mfilename
%   outputLocation
%
% variable input arguments, in order, character strings that are added to
% the diary file name (useful for tracking which diary file goes with
% which output file

[~,h] = system('hostname');
h = strtrim(h);
matD = now;
iD = datenum2iso(floor(matD),8); 
[~,~,~,hh,mm,ss] = datevec(matD);
ss = round(ss); %#ok<NASGU>
rng('shuffle')
D = [num2str(iD) 'T' num2str(hh,'%02d') num2str(mm,'%02d') '_' num2str(randi(intmax('int16')))];
for k=1:length(varargin)
    D = cat(2,D,'_',varargin{k});
end

output_file = [inputFunction '_' D '_' h '.txt'];

output_host = h;
if isdeployed
    diary(fullfile(outputLocation,output_file));
    disp(['host ' h ' starting Azure diary'])
    disp(['diary file ' fullfile(outputLocation,output_file)])
end
end