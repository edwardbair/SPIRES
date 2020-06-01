% WY=2012;
WY=2019;

vars={'fsca','fshade','grainradius','dust','weights','sensorZ'};
divisor=[100 100 1 10 100 1];
dtype={'uint8','uint8','uint16','uint16','uint8','uint8'};

for i=1:length(WY)
   matdates=datenum([WY(i)-1 10 1]):datenum([WY(i) 9 30]);
%  matdates=datenum([WY(i) 6 1]):datenum([WY(i) 6 20]);
%     [~,~,vars,divisor,dtype]=fill_and_run_modis(tiles,matdates,...
%         hdfbasedir,topodir,topofile,mask,R0,Ffile,shade,grain_thresh,dust_thresh,...
%         dustmask,tolval,outloc,nameprefix);

    out=smoothSPIREScube(nameprefix,vars,divisor,dtype,outloc,matdates,...
    nPersistDry,nPersistSnow,mingrainradius,maxgrainradius,mindust,maxdust,...
    mask,topofile,el_cutoff,fsca_thresh,cc,fice,...
    endconditions);

end

d=dir(fullfile(outloc,'*.h5'));
list=cell(length(WY),1);

for i=1:length(d)
    nm=d(i).name;
    for j=1:length(WY)
        wys=num2str(WY(j));
        rout=regexp(nm,['.*' wys '\.h5']);
        if rout
            list{j}=fullfile(outloc,nm);
        end
    end
end

make_modis_video(list,target,HUCunion,topofile);