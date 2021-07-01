function R=scaleMultiBandCube(R)
%scale
%scale multi-band cube to 1 as max value only if a band exceeds 1
%input sr mxnxb cube;
%output: scaled scube
mx=max(R,[],3);
t=any(R>1,3);
mx(~t)=1;
R=R./mx;