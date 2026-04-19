function PopObj = PF16OfZCAT(PopDec,K)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        % 也是一维流形，并且还是不连续的(K: 控制周期)。
    [N,M]=size(PopDec);
    M=M+1;

    PopObj=zeros(N,M);
    % 一维的流形
%     PopObj(:,1)=sin(PopDec(:,1)*pi/2);
% %     PopObj(:,2:M-2)=sin(PopDec(:,1)*pi/2).^(1+((2:M-2)-1)/(M-2));
%     PopObj(:,2:M-1)=sin(PopDec(:,1)*pi/2).^(1+((2:M-1)-1)/(M-2));
%     if(M>2)
%         PopObj(:,M-1)=1/2*(1+sin(5*pi*PopDec(:,1)-pi/2));      
%     end
%     PopObj(:,M)=(cos((2*K-1)*pi*PopDec(:,1))+2*PopDec(:,1)+4*K*(1-PopDec(:,1))-1)/(4*K);
% 
%     PopObj=(1:M).^(2).*PopObj;
%     % the sampled points 3-objective ZCAT  as follow:
%     % x=0:0.01:1
%     % [X,Y]=meshgrid(x);
%     % opt=PF1OfZCAT([X(:),Y
     PopObj(:,1)= sin(PopDec(:,1) *pi / 2);
   

    for j=2: M - 1
        PopObj(:,j) = sin(PopDec(:,1) *pi / 2).^(1.0 + (j - 1) / (M - 2));
%         assert 0 <= F[j - 1] <= 1.0
    end
    if (M > 2)
        PopObj(:,M-1) = 0.5 *(1 + sin(10 * PopDec(:,1) * pi / 2 - pi / 2));
    end
    PopObj(:,M) = (cos((2 * K - 1) * PopDec(:,1) *pi) + 2 * PopDec(:,1) + 4 * K * (1 - PopDec(:,1)) - 1) / (4 * K);
    PopObj=(1:M).^(2).*PopObj;
    
    
    
end

