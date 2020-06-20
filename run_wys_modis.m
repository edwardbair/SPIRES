WY=2017;
% WY=2001:2019;

for i=1:length(WY)
%    matdates=datenum([WY(i)-1 10 1]):datenum([WY(i) 9 30]);
  matdates=datenum([WY(i) 8 1]):datenum([WY(i) 8 30]);
%  fill_and_run_modis(tiles,matdates,...
%         hdfbasedir,topodir,topofile,mask,R0,Ffile,shade,grain_thresh,dust_thresh,...
%         dustmask,tolval,outloc,nameprefix);

    out=smoothSPIREScube(nameprefix,outloc,matdates,...
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

make_spires_video(list,target,HUCunion,topofile);