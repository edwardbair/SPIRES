function [cc,status]=create_cc(gfcc_dir,path,row,hdr)
%creates canopy cover image matching hdr for L8
%input:gfcc_dir - location of GFCC zip files, char
%path - path value, int
%row - row value, int
%hdr - hdr to mosaic and reproject to
%output: cc - canopy cover,0-1
%status - false, ok; true, error

gfcclist=cell(3,3);
pathoffset=-1:1;
rowoffset=-1:1;
status=false(size(gfcclist));

for i=1:length(pathoffset)
    for j=1:length(rowoffset)
        gfccname=fullfile(gfcc_dir,sprintf('GFCC30TC_p%03.0fr%03.0f_TC_2015.zip',...
            path+pathoffset(i),row+rowoffset(j)));
        try
            gfccfiles=unzip(gfccname,fullfile(gfcc_dir,'temp','gfcc'));
            gfcclist{i,j}=gfccfiles{1};
        catch
            fprintf('error w/ gfcc path:%i row:%i\n',path+pathoffset(i),row+rowoffset(j));
            status(i,j)=true;
        end
    end
end
gfcclist=gfcclist(:);
t=cellfun(@isempty,gfcclist);
gfcclist(t)=[];
    if ~isempty(gfcclist)
    [Big, BigRmap]=mosaicTiles(gfcclist,'tif','uint8',intmax('uint8'));
    %Big will come out as geog
    cc=rasterReprojection(Big,BigRmap,[],...
        hdr.ProjectionStructure,'rasterref',hdr.RasterReference,'method','nearest');
    cc=single(cc);
    cc(cc>100)=0;
    cc=cc./100;
    cc(cc<0.01)=0;
    else
        status=1;
        cc=[];
    end
end