function [ start,finish ] = contiguous(x)
% returns vectors (of equal length) start and finish, 
% correspoing to the start and finish of contiguously 
% positive blocks of x
% input: 
% x - vector of values
% output:
% start - N start indices of contig. blocks
% finsh - N finish indices of contig. blocks.
ind=find(x>0);
ind2=find(diff(ind)>1)+1;
start=ind(1);
start=[start; ind(ind2)];
finish=ind(end);
ind3=find(diff(x>0)<0);
finish=unique([ind3;finish]);

% t=finish-circshift(start,1);
% ind4=find(t<persist);
% if ~isempty(ind4);
%     ind4(ind4==1)=[];
%     ind4(ind4==length(finish))=[];
%     start(ind4)=[];
%     finish(ind4)=[];
% end

