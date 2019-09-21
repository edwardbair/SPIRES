
cm=cloudmask;
r=573;c=1522;
i=datenum([2010 4 10])-matdates(1)+1;
% i=66;
dt=datetime(matdates,'ConvertFrom','datenum');
figure;
hold on;
xx=squeeze(refl(r,c,:,:));

for j=1:7
    plot(dt,xx(j,:));
end


plot(dt,squeeze(cm(r,c,:)),'o');
plot(dt,squeeze(snowmask(r,c,:)),'x');
legend('b1','b2','b3','b4','b5','b6','b7','cloud','snow');

x=refl(:,:,:,i);
x(repmat(squeeze(cm(:,:,i)),[1 1 7]))=.3;
figure;image(x(:,:,[1 4 3]));axis image;
figure;image(refl(:,:,[1 4 3],i));axis image;