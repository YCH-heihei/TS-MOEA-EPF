function PopObj = PF14OfZCAT(PopDec)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        [N,M]=size(PopDec);
        M=M+1;
        
        PopObj=zeros(N,M);
        %% 一维的流形
        PopObj(:,1)=sin(pi/2*PopDec(:,1)).^2;
        PopObj(:,2:M-2)=sin(pi/2*PopDec(:,1)).^(2+((2:M-2)-1)/(M-2));
        if(M>2)
            PopObj(:,M-1)=1/2*(1+sin(3*pi*PopDec(:,1)-pi/2));
        end
        PopObj(:,M)=cos(pi/2*PopDec(:,1));
        
        PopObj=(1:M).^(2).*PopObj;
        % the sampled points 3-objective ZCAT  as follow:
        % x=0:0.01:1
        % [X,Y]=meshgrid(x);
        % opt=PF1OfZCAT([X(:),Y(:)]);
end

