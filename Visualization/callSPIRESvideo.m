
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

make_spires_video(list,target,HUCunion);