function MTL_list=MTL_parser(MTL_filename)
% MTL parser for LANDSAT datasets (MSS, TM, ETM+)
%
% Written in Nov 2012 by Evan Miles, Scott Polar Research Institute,
% University of Cambridge
%
% Extended from ETM+ MTL parser written by Seongsu Jeong,
% Geomatics and Remote Sensing Laboratory (GRSLAB), Yonsei University
%
% MAIN FUNCTIONALITY:
% Reads the information from MTL text file, which is provided with LANDSAT data.
%
% USAGE:
% MTL=MTL_parser(filename) - for a single MTL file
% MTL_list=MTL_parser(filename_list) - for multiple MTL file
%                                      The list should be in a matrix form
%                                      which each of its row vector is a
%                                      string of a MTL filename
% MTL_list=MTL_parser() - Searches all MTL files in current directory and
%                         parses them all
% 
% NOTIFICATION FOR THE OUTPUT VARIABLE:
% The values which surrounded "double quotation mark" was treated as a
% string, along with the information related to the time and date. Other
% variables were cosidered as a double-precision floating point variables.
%
% UPDATED FEATURES:
% 15-NOV-2012: Implemented recursive structure generator to interpret
%              nested groups with arbitrary names, removing necessity for 
%              identifying MTL-specific names.
% 08-MAR-2011: Capable of dealing with void input
%              When this case occurs, this module searches all MTL files
%              based on its filename and list them up. After then, this
%              function reads all detedted MTL files and returns the result
%              as a structure array in column vector form.
%
% LAST MODIFICATION: 15-NOV-2012
% 
%
%UPDATES SINCE DOWNLOAD:
%
%eval changed to evalc to suppress echoing of eval
%   supressed display in line 108


nMTL=0;
if nargin()==0 %when there is no input argument    
    file_list=ls;
    for cnt=1:size(file_list,1) %search all MTL files based on the file name
        if ~isempty(findstr(file_list(cnt,:),'MTL.txt'))
            nMTL=nMTL+1;
            MTL_filename(nMTL,:)=file_list(cnt,:);
        end
    end
elseif nargin()==1
    nMTL=size(MTL_filename,1);
else
    error('Incorrect input argument: Input a name of (or a list of) filename or make the input argument to void');
end

for cnt=1:nMTL
    % File open and line string input
    fin=fopen(MTL_filename(cnt,:),'r');
    str_in=fgetl(fin);

    %initialize structure parts
    group = {''};
    i=1;
    
    while ~strcmp(str_in,'END')
        % String parsing and refinement
        % input line refinement and processing
        
        [field value]=strtok(str_in,'=');
        
        %field name refinement
        %remove unnecessary space character from the field description
        field = regexprep(field,' ','');
        
        %Value refinement
        value=strtok(value,'=');
        
        %remove unnecessary space character from the field description
        while value(1)==' '
            value = value(2:end);
        end
        while value(end)==' '
            value = value(1:end-1);
        end
        
        %build structure
        if strcmpi(field,'GROUP')
            group{i} = value;
            i=i+1;
        elseif strcmpi(field,'END_GROUP')
            group = group(1:i-1);
            i=i-1;
        else
            if ~isempty(strfind(value,'"')) %if value contains ", replace with ' for matlab to eval as string
                value = regexprep(value,'"','''');
    %         if sum(value(1)=='"' && value(size(value,2))=='"') %If the value is a string wrapped by large quotation mark (")
    %             value=value(1,2:size(value,2)-1);
            elseif isempty(findstr(field,'TIME')) && isempty(findstr(field,'DATE')) %case for numeric data
                value=value;%no need to convert to number - should be as string for eval()
            else %for time or date data, wrap with '
                value = ['''' value ''''];
            end
            %disp(strcat(field,': ',value));
            j = 1;
            leadstring = '';
            while j<=length(group)
                leadstring = [leadstring group{j} '.'];
                j=j+1;
            end
            evalc([leadstring field '=' value]);
        end
        
        str_in=fgetl(fin);
    end
    
    evalc(['MTL_list(cnt,1)=' group{1} ';']);
    fclose(fin);
    
end
