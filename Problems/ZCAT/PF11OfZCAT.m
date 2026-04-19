function PopObj = PF11OfZCAT(PopDec,K)
%UNTITLED17 此处显示有关此函数的摘要
%   此处显示详细说明
        [N,M]=size(PopDec);
        M=M+1;
        
        PopObj=zeros(N,M);
        tempOfPopObj=cumsum(PopDec,2);
        PopObj(:,1)=1/(M-1)*tempOfPopObj(:,end);
        PopObj(:,2:end-1)=1./(M-(2:M-1)+1).*fliplr(tempOfPopObj(:,1:end-1)+(1-PopDec(:,end:-1:2)));
        PopObj(:,end)=(cos((2*K-1)*PopDec(:,1)*pi)+2*PopDec(:,1)+4*K*(1-PopDec(:,1)))/(4*K);
        
        PopObj=(1:M).^(2).*PopObj;
        % 第一个目标是前面所有变量的总和，其变量值应越小越好。
        % 而第二到倒数第二个目标则是希望在总和越小的情况下，变量值越大越好。
        % 而最后一个目标是希望第一个变量越大越好，但是这里cos((2K-1)*PopDec(:,1))则是起到了让这条1-PopDec的直线存在+1的波动。
        % 实际上的波动本身就是造成不连续的原因，因为直线也仅仅只是一条直线由1到0的直线，直线本身就是一个PF。
        % 直线的斜率是
        % 减少的比增加的多是造成不连续的原因。
        % cos的增加值随x的周期增加有时会比直线中y随x的减少值要大，这就出现了x较小时对应的y比x较大时对应y更大的情况。
        % 不不不，这一点不能保证不连续。
        % cos可能是减小也可能是增加，在减小时将直线拉低。
        % cos增加时又可能将直线拉高。
        % 这一低一高之间便是形成不连续PF的区域。
        % 对于cos的周期: 前pi都是有效区域，接着的后
        
        
        % the sampled points 3-objective ZCAT  as follow:
        % x=0:0.01:1
        % [X,Y]=meshgrid(x);
        % opt=PF1OfZCAT([X(:),Y(:)]);
end

