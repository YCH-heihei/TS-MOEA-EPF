function PopObj = PF8OfZCAT(PopDec)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        [N,M]=size(PopDec);
        M=M+1;
        
        PopObj=zeros(N,M);
        tempOfPopObj=cumprod(1-sin(PopDec.*pi/2),2);
        PopObj(:,1)=1-tempOfPopObj(:,end);
        PopObj(:,2:end-1)=1-fliplr(tempOfPopObj(:,1:end-1).*(1-cos(PopDec(:,end:-1:2)*pi/2)));
        PopObj(:,end)=cos(pi/2*PopDec(:,1));
        PopObj=(1:M).^(2).*PopObj;
        % the sampled points 3-objective ZCAT  as follow:
        % x=0:0.01:1
        % [X,Y]=meshgrid(x);
        % opt=PF1OfZCAT([X(:),Y(:)]);
end

