function PopObj = PF5OfZCAT(PopDec)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        [N,M]=size(PopDec);
        M=M+1;
        
        PopObj=zeros(N,M);
        PopObj(:,1:M-1)=PopDec(:,1:M-1);
        PopObj(:,M)=(exp(1/(M-1).*sum(1-PopDec,2)).^(8)-1)/(exp(1).^8-1);
        PopObj=(1:M).^(2).*PopObj;
        % the sampled points 3-objective ZCAT  as follow:
        % x=0:0.01:1
        % [X,Y]=meshgrid(x);
        % opt=PF1OfZCAT([X(:),Y(:)]);
end

