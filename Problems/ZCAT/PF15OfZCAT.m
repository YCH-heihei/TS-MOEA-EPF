function PopObj = PF15OfZCAT(PopDec,K)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        % 也是一维流形，并且还是不连续的(K: 控制周期)。
        [N,M]=size(PopDec);
        M=M+1;
        
        PopObj=zeros(N,M);
        %% 一维的流形
        PopObj(:,1:M-1)=PopDec(:,1).^(1+(1:M-1)/(4*M));
        PopObj(:,M)=(cos((2*K-1)*pi*PopDec(:,1))+2*PopDec(:,1)+4*K*(1-PopDec(:,1))-1)/(4*M);
        
        
        PopObj=(1:M).^(2).*PopObj;
        % the sampled points 3-objective ZCAT  as follow:
        % x=0:0.01:1
        % [X,Y]=meshgrid(x);
        % opt=PF1OfZCAT([X(:),Y(:)]);
end

