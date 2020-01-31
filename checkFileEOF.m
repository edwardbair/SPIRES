function checkFileEOF(filename)
% some weird code to make sure EOF for file, based on dump of
% input file that smarts295.exo or smarts295bat.exe would not open because
% of lack of EOF
%Read file, adds EOF if needed, and rewrites file

fileID = fopen(filename,'r','n');
if fileID==-1
    warning('unable to open file ''%s'', returning',filename)
else
    X = uint8(fread(fileID,'uint8'));
    fclose(fileID);
    newline = uint8(10);
    eomark = uint8(13);
    k = find(X==newline);
    if ~isempty(k)
        k = flipud(k)';
        % check to make sure that we need to do this
        if X(k(1)-1)~=eomark
            % replace all newline with [eomark newline]
            for n=k
                newX = cat(1,X(1:n-1),eomark,X(n:end));
                X=newX;
            end
            % rewrite the file, in binary
            delete(filename);
            fileID = fopen(filename,'w');
            fwrite(fileID,X);
            fclose(fileID);
        end
    end
end
end
