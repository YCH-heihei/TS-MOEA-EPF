classdef ZCAT9 < PROBLEM
% <multi/many> <real> <large/none> <expensive/none>
% Benchmark MOP proposed by Walking Fish Group
% K --- 3 --- The position parameter, which should be a multiple of M-1
% Bias --- 0 --- bias
% imbalance --- 1 --- balance
% complicated --- 0 ---  Complicated PS flag
% level --- 1 ---  Complicated PS flag

%------------------------------- Reference --------------------------------
% S. Huband, P. Hingston, L. Barone, and L. While, A review of
% multiobjective test problems and a scalable test problem toolkit, IEEE
% Transactions on Evolutionary Computation, 2006, 10(5): 477-506.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(Access = private)
        K;  % Position parameter.
        m; % num of Position-related decision variables;
        Bias; 
        imbalance; % different distance functions for different objectives;
        level; % difficulty approaching the true PF, the higher the difficulty, the more difficult it is;
        complicatedParetoSet; % The method of associating distance variables with position variables;
    end
    methods
        %% Default settings of the problem
        function Setting(obj)
            if isempty(obj.M); obj.M = 3; end
            obj.m=obj.M-1;
            [obj.K,obj.Bias,obj.imbalance,obj.complicatedParetoSet,obj.level] = obj.ParameterSet(3,0,1,0,1);
%             if isempty(obj.D); obj.D = obj.M*10; end
%             obj.m = obj.ParameterSet(obj.M-1);% 决策变量的维度数。
%             obj.Bias = obj.ParameterSet();
%             obj.imbalance = obj.ParameterSet(0);
%             obj.level = obj.ParameterSet(1);
%             obj.complicatedParetoSet = obj.ParameterSet(0);
%             obj.K=obj.ParameterSet(3);
            if isempty(obj.D); obj.D = obj.M*10; end
            obj.lower    = -(1:obj.D)/2;
            obj.upper    = (1:obj.D)/2;
            obj.encoding = ones(1,obj.D);
        end
        %% Calculate objective values
        function PopObj = CalObj(obj,PopDec)
            % 这里是首先标准化的，这使得整个搜索空间依然是在
            PopDec=(PopDec-obj.lower)./(obj.upper-obj.lower);
            m=obj.m;
            n=obj.D;
            M=obj.M;
            level=obj.level;
            % 前M个目标函数。
            
            PopObj=PF9OfZCAT(PopDec(:,1:m));
            
%             PopObj=PopObj.*((1:M).^2);
            
            % 对收敛变量借助g函数做偏移。
            if(obj.complicatedParetoSet==0)
                PopDec(:,m+1:n)=PopDec(:,m+1:n)-g0(PopDec(:,1:m),n);
            else
                PopDec(:,m+1:n)=PopDec(:,m+1:n)-g5(PopDec(:,1:m),n);
            end
            
            if(obj.Bias==1)
                PopDec(:,m+1:n)=ZBias(PopDec(:,m+1:n));
            end
            % 对变量基于模M将其分配给各个变量上。
            indexOfM=mod((m+1:n)-m,obj.M);
            
            indexOfM(indexOfM==0)=obj.M;
            for i=1:obj.M
                if(obj.imbalance)
                    if(mod(i,2)==0)
%                         PopObj(:,i)=PopObj(:,i)+i^2*Z4(PopDec(:,m+find(indexOfM==i)),obj.K);
                        PopObj(:,i)=PopObj(:,i)+i^2*Z3(PopDec(:,m+find(indexOfM==i)),5);
                    else
                        PopObj(:,i)=PopObj(:,i)+i^2*Z1(PopDec(:,m+find(indexOfM==i)),5);
                    end
                else
                    temp=eval(strcat('Z',num2str(obj.level),'(PopDec(:,m+find(indexOfM==i)),obj.K);'));
                    PopObj(:,i)=PopObj(:,i)+i.^2*temp;
                end
            end
        end
        
        %% Generate points on the Pareto front
        function R = GetOptimum(obj,N)
         % the sampled points 3-objective ZCAT  as follow:
            if(obj.M==3||obj.M==2)
                datePath='.\Problems\Multi-objective optimization\ZCAT';
                tempPath=strcat(datePath,'\','ZCAT9','\','ZCAT9','_M',num2str(obj.M),'_Num',num2str(N),'.mat');
                load(tempPath);
                R=opt;
%             elseif(obj.M==2)
%                 x=0:0.01:1;
%                 opt=PF9OfZCAT(x');
%                 [Fno,~]=NDSort(opt,1);
%                 opt=opt(Fno==1,:);
%                 opt=sortrows(opt,1);
%                 R=opt;
            else
                R=[];
            end
            
            % 从4万个点种反复删除
        end
        %% Generate the image of Pareto front
        function R = GetPF(obj)
            if(obj.M==3)
%                 datePath='.\Problems\Multi-objective optimization\ZCAT';
%                 tempPath=strcat(datePath,'\','ZCAT1','_M',num2str(obj.M),'_Num',num2str(441),'.mat');
%                 load(tempPath);
                x=0:0.05:1;
                [X,Y]=meshgrid(x);
                opt=PF9OfZCAT([X(:),Y(:)]);
                [Fno,~]=NDSort(opt,1);
                opt=opt(Fno==1,:);
%                 firstAxis=zeros(1,size(opt,2));
%                 firstAxis(1)=1;
%                 zMin=min(opt,[],1);
%                 zMax=max(opt,[],1);
%                 normOfOpt=(opt-zMin)./(zMax-zMin);
%                 direction=normOfOpt./sum(normOfOpt,2);
%                 cosineDis=pdist2(direction,firstAxis,'cosine');
%                 [~,orderIndex]=sort(cosineDis);
%                 opt=opt(orderIndex,:);
                num=floor(size(opt,1).^(1/2));
                R={reshape(opt(:,1),num,[]),reshape(opt(:,2),num,[]),reshape(opt(:,3),num,[])};
            elseif(obj.M==2)
                x=0:0.01:1;
                opt=PF9OfZCAT(x');
                [Fno,~]=NDSort(opt,1);
                
                opt(Fno~=1,:)=nan;
                opt=sortrows(opt,1);
                R=opt;
            else
                R=[];
            end
        end
    end
end
