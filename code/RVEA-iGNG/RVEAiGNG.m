classdef RVEAiGNG < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% RVEA based on improved growing neural gas
% alpha --- 2 --- The parameter controlling the rate of change of penalty

%------------------------------- Reference --------------------------------
% Q. Liu, Y. Jin, M. Heiderich, T. Rodemann, and G. Yu, An adaptive
% reference vector-guided evolutionary algorithm using growing neural gas
% for many-objective optimization of irregular problems, IEEE Transactions
% on Cybernetics, 2022, 52(5): 2698-2711.
%--------------------------------------------------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

% This function is written by Qiqi Liu

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            params.N = Problem.N;% 种群规模.
            params.MaxIt = 50;% GNG训练次数。
            params.L = 50;  % 产生新节点的频率。
            params.epsilon_b = 0.2; % 最优节点的学习率。
            params.epsilon_n = 0.006; % 邻居节点的学习率。
            params.alpha = 0.5; % 产生子代所使用的折扣。
            params.delta = 0.995; % 节点累积误差的衰减。
            params.T = 50; % age未更新的上界。
            
            %% Parameter setting
            alpha = Algorithm.ParameterSet(2);
            [V,Problem.N] = UniformPoint(Problem.N,Problem.M);% 初始V为均匀的参考向量。
            Population    = Problem.Initialization();
            net = InitilizeGrowingGasNet(V,Population,params);% 使用当前种群及有效参考向量来对网络进行初始化。
            Archive = UpdateArchive(Population,[],Problem.N);
            scale = ones(1,Problem.M); % 先让scale为。
            zmin = min(Population.objs,[],1); % 理想点。
            genFlag = []; % 

            while Algorithm.NotTerminated(Population)
                MatingPool = randi(length(Population),1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                zmin       = min([zmin;Offspring.objs],[],1); % 更新理想点。
                % 对于MaF6问题，RVEAiGNG甚至刚开始仅仅选择几个个体,这是因为RVEAiGNG仅仅从关联的参考向量中选择个体。
                [Population,net,V,Archive,scale,genFlag] = EnvironmentalSelection([Population,Offspring],V,(Problem.FE/Problem.maxFE)^alpha,net,params,Archive,Problem,scale,zmin,genFlag);
            end
            temp=0;
        end
    end
end