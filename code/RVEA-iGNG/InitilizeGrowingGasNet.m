function net = InitilizeGrowingGasNet(V,Population,params)

%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    N = params.N;
    MaxIt = params.MaxIt;% 最大迭代次数。
    L = params.L; % 产生新节点的频率。
    epsilon_b = params.epsilon_b;% 最近节点学习率。
    epsilon_n = params.epsilon_n;% 邻居节点学习率。
    alpha = params.alpha; % 产生新节点的折扣。
    delta = params.delta; % 其它节点累积误差的衰减率。
    T = params.T; % 边的最大邻居。

    PopObj = Population.objs;
    [NP,M]  = size(PopObj);
    PopObj = PopObj - repmat(min(PopObj,[],1),NP,1);% 这里只是原点化，没有进行标准化。
    Angle = acos(1-pdist2(PopObj,V,'cosine')); % 这里使用夹角进行关联。
    [~,associate] = min(Angle,[],2);
    valid = unique(associate); % 使用夹角来将个体关联到初始的参考向量上。
    RefSize = size(valid,1);% 有效参考向量的个数。
    %% Initialization
    Ni = 2;
    w = zeros(Ni, M);
    if RefSize>=2
        % 先取前两个有效的参考向量。
        for i = 1:Ni
            w(i,:) = V(valid(i,:),:);
        end
    else
        w(1:Ni,:) = V(randperm(N,Ni),:);% 若没有无效的，则随机选择两个。
    end
    E = zeros(Ni,1);
    C = zeros(Ni, Ni);
    t = zeros(Ni, Ni);
    nx = 0;
    %% Loop
    for it = 1:MaxIt
        for kk = 3:RefSize
            % Select Input
            nx = nx + 1;
            x = V(valid(kk,:),:); % 有效参考向量中的第kk个。

            % Competion and Ranking
            d = pdist2(x, w);
            [~, SortOrder] = sort(d);
            s1 = SortOrder(1); % 最近个体。
            s2 = SortOrder(2); % 次近个体。

            % Aging(更新边age,表示存在数据关联，并用于删除旧边)。
            t(s1, :) = t(s1, :) + 1; 
            t(:, s1) = t(:, s1) + 1;

            % Add Error
            E(s1) = E(s1) + d(s1)^2;

            % Adaptation
            w(s1,:) = w(s1,:) + epsilon_b*(x-w(s1,:)); % 以epsilon_b更新最近节点。
            Ns1 = find(C(s1,:)==1); % 最近节点邻居。
            % 以epsilon_nb更新最近节点的邻居。
            for j=Ns1
                w(j,:) = w(j,:) + epsilon_n*(x-w(j,:));
            end

            % Create Link
            C(s1,s2) = 1;
            C(s2,s1) = 1;
            t(s1,s2) = 0;
            t(s2,s1) = 0;

            % Remove Old Links
            C(t>T) = 0;
            nNeighbor = sum(C);
            % 删除孤立节点。
            AloneNodes = (nNeighbor==0);
            C(AloneNodes, :) = [];
            C(:, AloneNodes) = [];
            t(AloneNodes, :) = [];
            t(:, AloneNodes) = [];
            w(AloneNodes, :) = [];
            E(AloneNodes) = [];

            % Add New Nodes
            if mod(nx, L) == 0 && size(w,1) < N
                [~, q] = max(E);
                [~, f] = max(C(:,q).*E);
                r = size(w,1) + 1;
                w(r,:) = (w(q,:) + w(f,:))/2;
                C(q,f) = 0;
                C(f,q) = 0;
                C(q,r) = 1;
                C(r,q) = 1;
                C(r,f) = 1;
                C(f,r) = 1;
                t(r,:) = 0;
                t(:, r) = 0;
                E(q) = alpha*E(q);
                E(f) = alpha*E(f);
                E(r) = E(q);
            end

            % Decrease Errors
            E = delta*E;
        end
    %     PlotResults(w, C)
    end

    for ii = 1:size(w,1)
        ageSum(ii,:) = sum(t(ii,find(C(ii,:) == 1),:),2);
        ageSumBefore = ageSum; % 初始化ageSumBeform，用以更新flag。
        flag(ii,:) = 0;
    end
    net.w = w;
    net.E = E;
    net.C = C;
    net.t = t;
    net.nx = nx;
    net.ageSumBefore = ageSumBefore;
    net.flag = flag;
end