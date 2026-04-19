function PopObj = PF6OfZCAT(PopDec)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        [N,M]=size(PopDec);
        M=M+1;
        
        k=40;
        rou=0.05;
        miu=1/(M-1).*(sum(PopDec,2));
        
        PopObj=zeros(N,M);
        PopObj(:,1:M-1)=PopDec(:,1:M-1);
        PopObj(:,M)=((1+exp(2*k*miu-k)).^(-1)-rou*miu-(1+exp(k)).^(-1)+rou)...
            /((1+exp(-k)).^(-1)-(1+exp(k)).^(-1)+rou);
        PopObj=(1:M).^(2).*PopObj;
        % the sampled points 3-objective ZCAT  as follow:
        % x=0:0.01:1
        % [X,Y]=meshgrid(x);
        % opt=PF1OfZCAT([X(:),Y(:)]);
end

