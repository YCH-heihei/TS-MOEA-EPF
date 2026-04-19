classdef RVEAAdaK < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation> <constrained/none>
% Reference vector guided evolutionary algorithm
% alpha ---   2 --- The parameter controlling the rate of change of penalty
% fr    --- 0.1 --- The frequency of employing reference vector adaptation

%------------------------------- Reference --------------------------------
% R. Cheng, Y. Jin, M. Olhofer, and B. Sendhoff, A reference vector guided
% evolutionary algorithm for many-objective optimization, IEEE Transactions
% on Evolutionary Computation, 2016, 20(5): 773-791.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [alpha,fr] = Algorithm.ParameterSet(2,0.1);

            %% Generate the reference points and random population
            [W,Problem.N] = UniformPoint(Problem.N,Problem.M);
            Population     = Problem.Initialization();
            % V              = V0;

            V=W;
            fstep=0.05;
            flast=0.9;
            A=Population.objs; % 源代码及原文中都是直接选择
            Wadapt=W;
            
            zmin=min(A,[],1);
            flag=6; % 1: RSE, 2: GAE, 3: COU, 4:PT, 5:MPT, 6: KRA，从图中得到RVEA使用KRA较优。
            znadir=[]; % 极差点。
            scale=false; % 是否需对参考向量依据目标值进行缩放(RVEA是在外侧进行的缩放)。
            update=false; % 是否需要遵守第i个个体在第i个参考向量上表现最好(RVEA因未使用小生境，从而不需要采取基于相邻个体来产生配对池)。
            T=[]; % 小生境个数(对于RVEA则不需要)。
            B=[]; % 小生境下标(对于RVEA则不需要)。
            utility=[]; % 参考向量上的聚合函数。
            
            % zmax=max(A(NDSort(A,1)==1,:),[],1);
            r=ones(1,Problem.M);
            % r(r==0)=1;
            %% Optimization
            while Algorithm.NotTerminated(Population)
                MatingPool = randi(length(Population),1,Problem.N);
                Offspring  = OperatorGA(Problem,Population(MatingPool));
                zmin=min(zmin,min(Offspring.objs,[],1));
                Population = EnvironmentalSelectionRVEAAdaK([Population,Offspring],zmin,V,(Problem.FE/Problem.maxFE)^alpha);
                % if ~mod(ceil(Problem.FE/Problem.N),ceil(fr*Problem.maxFE/Problem.N))
                %     V(1:Problem.N,:) = ReferenceVectorAdaptation(Population.objs,V0);
                % end
               % if ~mod(ceil(Problem.FE/Problem.N),ceil(fr*Problem.maxFE/Problem.N))
               %      V(1:Problem.N,:) = ReferenceVectorAdaptation(Population.objs,V0);
               %  end
               if mod(ceil(Problem.FE/Problem.N), ceil(fr*Problem.maxFE/Problem.N)) == 0
                    r = referenceVectorAdaptation(Population, zmin);
               end
                          
               [A, Wadapt,~,~] = AdaK(A, Offspring, Population, Wadapt, W, Problem.N, zmin, fstep, flast, ceil(Problem.FE/Problem.N), ceil(Problem.maxFE/Problem.N), flag, znadir, scale, update, T, B, utility);
               V = Wadapt.*r;
            end
        end
    end
end

