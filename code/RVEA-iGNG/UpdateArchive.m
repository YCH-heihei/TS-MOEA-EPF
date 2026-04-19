function Archive = UpdateArchive(Population,Archive,MaxSize)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    Archive = [Archive,Population];
    o = Archive.objs;
    [c,ia,ic] = unique(o,'rows');
    Archive = Archive(ia);
    ND = NDSort(Archive.objs,1);
    Archive = Archive(ND==1);
    N  = length(Archive);
    if N <= MaxSize
        return;
    end
    
    %% Calculate the fitness of each solution
    CAObj = Archive.objs;
    % 使用当前的最小值及最大值来进行初始化。
    CAObj = (CAObj-repmat(min(CAObj),N,1))./(repmat(max(CAObj)-min(CAObj),N,1));
    I = zeros(N);
    for i = 1 : N
        for j = 1 : N
            I(i,j) = max(CAObj(i,:)-CAObj(j,:)); % 个体j相对于个体i的最大优势。
        end
    end
    C = max(abs(I)); 
    %  个体j相对于其它个体的最大优势中的最大优势，如果节点存在某个值最小， 那么这个值就是1。
    % 如果节点所有目标的中心，那么这个值接近0.5。
    F = sum(-exp(-I./repmat(C,N,1)/0.05)) + 1; % 最后的exp相当于累乘。
    
    %% Delete part of the solutions by their fitnesses
    Choose = 1 : N;
    while length(Choose) > MaxSize
        [~,x] = min(F(Choose)); % 反复删除F值最小的个体。
        F = F + exp(-I(Choose(x),:)/C(Choose(x))/0.05); % 重新计算F值。
        Choose(x) = [];
    end
    Archive = Archive(Choose);

    %% Delete those which are too far from the archive
    o = Archive.objs;
    o = o-repmat(min(o),size(o,1),1);
    d = sqrt(sum(o.^2,2));
    meanD = sum(d,1)/size(o,1);
    delete = find(d>10*meanD);
    Archive(delete) = [];
    % 注意，它的存档并没有强制保留极大值点。
end