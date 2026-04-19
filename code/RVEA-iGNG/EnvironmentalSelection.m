function [Population,net,V,Archive,scale,genFlag] = EnvironmentalSelection(Population,V,theta,net,params,Archive,Problem,scale,zmin,genFlag)
% The environmental selection of RVEA

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu
    % Population: 包含当前子代及种群。
    % theta: 当前的进度。
    % params: 网络有关的参数。
    % 当前更新的scale值。
    Population = Population(NDSort(Population.objs,1)==1); % 仅仅保留非被支配个体。
    PopObj = Population.objs;
    [N,M]  = size(PopObj);
    %% Translate the population
    PopObj = PopObj - repmat(zmin,N,1); % 原点化。
    %% delete the outliers
    d = sqrt(sum(PopObj.^2,2)); % 到原点的距离。
    meanD = sum(d,1)/size(PopObj,1);% 到原点距离的均值。 
    delete = find(d>10*meanD); % 认为这些点是离群点(这一点可以考虑)
    Population(delete)=[];
    PopObj(delete,:)=[]; % 删除离群点。
    SavePopObj = Population.objs;% 当前子代加入后的非被支配个体。


    Archive = UpdateArchive(Population,Archive,2*Problem.N); % 种群规模为2N。
    ArcObj = Archive.objs;

    wholeObj = [SavePopObj;ArcObj];% 整个非被支配个体，即使未被参与到种群当中，也将使用其用于训练，但是。
    Population = [Population Archive]; % 注意这里包含了决策变量，而不仅仅是目标空间，所以可以用来引导种群的筛选。
    % 因此便于根据需要取回某个目标向量对应的决策变量。
    [c,ia,ic] = unique(wholeObj,'rows');
    Population = Population(ia);
    wholeObj = wholeObj(ia,:);
    wholeObj = wholeObj - repmat(zmin,size(wholeObj,1),1); % 原点化。



    wholeObj1 = wholeObj./scale;% 标准化，scale: 采取定期更新的方法。
    temp1 = wholeObj1./sum(wholeObj1,2); % 平面化。


    PopObj = wholeObj;
    [N,M]  = size(PopObj);
    fr = 0.1;

    gen    = ceil(Problem.FE/Problem.N);% 当前代数。
    maxgen = ceil(Problem.maxFE/Problem.N); % 最大代数。
    if ~mod(gen,ceil(fr*maxgen))&& gen <= round(1*maxgen)
        % 间隔0.1阶段将更新跨度并且处于可更新的范围内。
        scale = max(ArcObj,[],1)-min(ArcObj,[],1);% 更新跨度。
    end
    % temp1: 当前非被支配—位于种群、子代、存档当中的。
    if size(temp1,1) > 2&& isempty(find(isnan(temp1)==true))&& gen <= round(1*maxgen) && isempty(genFlag)
        % 使用投影在超平面上的存档来更新GNG网络。
        % 注意这里的V不在一定为初始的参考向量了，而是上一代中用以引导种群优化的参考向量。
        % temp1: 当前非被支配个体在超平面上的投影。
        % [[];ArcObj]: 当前筛选后的存档对应的目标向量。
        % genFlag: 初始为空，也仅仅在为空时，需要生成产考向量。
         % 这个算法也是每代就更新的，那为什么它能够摆脱极值点？
        [V,net,genFlag] = TrainGrowingGasNet(V,temp1,net,scale,params,Problem,[[];ArcObj],genFlag,zmin);
    end
    NV     = size(V,1); 

    %% Calculate the degree of violation of each solution

    CV = sum(max(0,Population.cons),2);%(违反约束程度)

    %% Calculate the smallest angle value between each vector and others
    cosine = 1 - pdist2(V,V,'cosine');
    cosine(logical(eye(length(cosine)))) = 0;
    gamma  = min(acos(cosine),[],2); % 参考向量与其它参考向量的最小夹角，用以作为APD值标准化所使用的值。

    %% Associate each solution to a reference vector
    Angle = acos(1-pdist2(PopObj,V,'cosine'));
    % 使用余弦值将个体关联到与其最近的参考向量上，注意这里的V仅仅是原点化后的，并未缩放。
    [~,associate] = min(Angle,[],2);

    %% Select one solution for each reference vector
    Next = zeros(1,NV); % NV: 参考向量数目，为每个参考向量关联其最优APD值的个体。
    for i = unique(associate)' % 针对对每个存在个体关联的参考向量。
        current1 = find(associate==i & CV==0);% 不违背约束的个体。
        current2 = find(associate==i & CV~=0);% 违背约束的个体。
        if ~isempty(current1)
            % Calculate the APD value of each
            % solution(夹角比上gamma值)乘上个体到原点的欧式距离(标准化后的)。
            APD = (1+M*theta*Angle(current1,i)/gamma(i)).*sqrt(sum(PopObj(current1,:).^2,2));
            % Select the one with the minimum APD value
            [~,best] = min(APD);
            Next(i)  = current1(best);
            % 为每个参考向量关联其APD值最好的个体。
        elseif ~isempty(current2)
            % Select the one with the minimum CV value
            [~,best] = min(CV(current2)); % 选择违背约束最小的个体。
            Next(i)  = current2(best);
        end
    end
    % Population for next generation
    Population1 = Population(Next(Next~=0)); % 选择之前已经选择的个体。
    %% select the corner solutions in each generation

    fm = [];
    selectedFirst = unique([Next(Next~=0) fm]);
    Population = [Population(selectedFirst)];

    if length(Population) > Problem.N 
        % 在种群规模大于N时进一步进行筛选(因为可能参考向量数目大于N)，但选择的个体也仅仅是从已选则个体中进行筛选。
        PopObj = Population.objs;
        zmax = max(PopObj,[],1);
        PopObj = PopObj - repmat(zmin,size(PopObj,1),1); % 仅仅是原点化。

        temp1 = PopObj;

        Choose = false(1,size(temp1,1));
        [~,Extreme1] = min(temp1,[],1);
        [~,Extreme2] = max(temp1,[],1);
        Choose(Extreme1) = true; % 优先保留极小值点。
        Choose(Extreme2) = true; % 优先保留极大值点。

        while sum(Choose) < Problem.N
            ind = find(Choose== false);% 未选择。
            choId = setdiff(1:size(temp1,1),ind); % 已经选择。
            PopObj1temp = temp1(choId,:);     
            WholeObjtemp = temp1(ind,:);
            dis =  pdist2(WholeObjtemp,PopObj1temp); % 未选择个体到已选择个体的距离。
            [mindis,] = min(dis,[],2);
            [~,associate] = max(mindis,[],1); % 选择与已经选择个体最不拥挤的个体。
            Choose(ind(associate)) = true;
        end
        Population = Population(Choose);
    end
end