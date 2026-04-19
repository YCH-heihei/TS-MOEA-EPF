function PopObj = PF17OfZCAT(PopDec)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        % 也是一维流形，并且还是不连续的(K: 控制周期)。
    [N,M]=size(PopDec);
    M=M+1;

    PopObj=zeros(N,M);
    %% 一维的流形
    index=all(PopDec<=0.5,2);
    PopObj(index,1:M-1)=repmat(PopDec(index,1),1,M-1);
    PopObj(~index,1:M-1)=PopDec(~index,1:M-1);
    PopObj(index,M)=(exp(1-PopDec(index,1)).^(8)-1)/(exp(1).^8-1);
    PopObj(~index,M)=(exp(1/(M-1)*sum(1-PopDec(~index,:),2)).^8-1)/(exp(1).^8-1);
    
    PopObj=(1:M).^(2).*PopObj;
    % the sampled points 3-objective ZCAT  as follow:
    % x=0:0.01:1
    % [X,Y]=meshgrid(x);
    % opt=PF1OfZCAT([X(:),Y(:)]);
end

