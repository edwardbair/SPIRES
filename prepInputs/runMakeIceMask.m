function runMakeIceMask(S,list,outloc)
%run makeIceMask on a list of files w/ hdr struct
%input: S - GLIMS shapefile (no M & no Z values)
%list list of 
for i=1:length(list)
    t1=tic;
    m=matfile(list{i});
    fice=makeIceMask(S,m.hdr,3);
    [~,name]=fileparts(m.Properties.Source);
    s=strsplit(name,'_');
    fname=fullfile(outloc,[s{1},'.mat']);
    mm=matfile(fname);
    mm.fice=fice;
    mm.hdr=m.hdr;
    t2=toc(t1);
    fprintf('wrote %s in %g min\n',fname,t2/60);
end