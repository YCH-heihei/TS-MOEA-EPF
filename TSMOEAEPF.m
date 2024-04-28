classdef TSMOEAEPF < ALGORITHM
% <multi/many> <real/integer/label/binary/permutation>
% Decomposition based evolutionary algorithm guided by growing neural gas
% aph ---   0.1  --- Parameter alpha
% eps --- 0.314 --- Parameter epsilon


%------------------------------- Reference --------------------------------
% Y. Liu, H. Ishibuchi, N. Masuyama, and Y. Nojima, Adapting reference
% vectors and scalarizing functions by growing neural gas to handle
% irregular Pareto fronts. IEEE Transactions on Evolutionary Computation,
% 2020, 24(3): 439-453.
%--------------------------------------------------------------------------
% Copyright Yiping Liu
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------

    methods
        function main(Algorithm,Problem)
            %% Parameter setting
            [aph,eps] = Algorithm.ParameterSet(0.1,0.314); 
            
            %% DEA Initialization
            [Ru,Problem.N] = UniformPoint(Problem.N,Problem.M);	% Uniform reference vectors
            Population     = Problem.Initialization();          % Random population
            zMin           = min(Population.objs,[],1);         % Ideal Point  
            AS             = [];                                % Input Signal Archive
            Ruq            = Ru;                                % Refernce vectors in Ru for selection
            PopObjs=Population.objs;
            [FrontNo,~]    = NDSort(Population.objs,Problem.N); % First Fitness for the first mating selection
            firstIndex=find(FrontNo==1);
            PopObjs=PopObjs(firstIndex,:);
            AS=Population(firstIndex);
            M=Problem.M;
%             AS=[PopObjs];
            zMax=max(PopObjs,[],1);
            crd            = zeros(1,Problem.N);                % second Fitness for the first mating selection   
            MaxGen         = ceil(Problem.maxFE/Problem.N);    	% Maximum Generation
%             FrontNo = NDSort(Population.objs,Problem.N);
            
            %% GNG Initialization
%             spand=floor(1/20*Problem.M*Problem.N);
%             spand=ceil(1/4*Problem.N);      
%             ArchiveSize1 = ceil(1/2*Problem.M)*Problem.N;
%             ArchiveSize2 = ceil(1/2*Problem.M)*Problem.N;
           ArchiveSize1 = ceil(2)*ceil(Problem.N+3);
            ArchiveSize2 = ceil(2)*ceil(Problem.N+3);
            NoG = ceil(aph*MaxGen);                   % Number of generations of Not Training GNGNumber of iterations to train GNG per Generation(针对每轮训练需要对数据进行多少次重复训练)  
            GNGnet.maxIter = 1;                 % Number of iterations to train GNG per Generation  
%             GNGnet.maxAge = (Problem.N);          % Maximum cluster age 
            GNGnet.maxAge = (10*Problem.M); 
%             GNGnet.maxAge = (Problem.M).^2;
            GNGnet.maxNode = ceil(1*Problem.N);         % Max number of nodes
            GNGnet.lambda = 0.05*Problem.N;      % Cycle for topology reconstruction   
            GNGnet.hp = [];                     % Hit point of node 
            GNGnet.maxHP = 2*ArchiveSize2;       % Max HP of node 
            GNGnet.Node = [];                   % Node
            GNGnet.NodeS = [];                  % Expanded node 
            GNGnet.NodeP = [];                  % Node mapped to hyperplane 
            GNGnet.Err = zeros(2,2);            % Error(修改，将误差累积到边，而不是节点上) 
            GNGnet.edge = zeros(2,2);           % Edge between nodes 
            GNGnet.age = zeros(2,2);            % Age of edge 
            GNGnet.epsilon_a = 0.2;             % Learning coefficient
            GNGnet.epsilon_nb = 0.02;           % Learning coefficient of neighbor
            GNGnet.alpha = 0.5;                 % Nodes r1max and r2max error reduction constant
            GNGnet.delta = 1/5.^(1/ArchiveSize1);                 % Error reduction coefficient  
            trainedFlag=0;
            EM=eye(Problem.M);
            %% Optimization
            while Algorithm.NotTerminated(Population)


              %% 原始方案。
                                          
              
%               % GNG-based adaptation(当前代属于可训练GNG的代)
                if ceil(Problem.FE/Problem.N) <= MaxGen - NoG
                    %% 原始方案。
                    MatingPool=TournamentSelection(2,Problem.N,FrontNo);
                    FirstIndex=find(FrontNo==1);
                    PopObjs=Population.objs;
                    objs=PopObjs(FirstIndex,:);
                    [zMax,maxIndex]=max(objs,[],1);


                    Offspring  = OperatorGA(Problem,Population(MatingPool));
                     if(Problem.FE>=16340)
                            % rng(20000);
                            temp=0;
                     end
                     if(Problem.FE>=209)
                            temp=0;
                     end
                    zMin= min([zMin;Offspring.objs],[],1);
                    % Input Signal Archive Update
                    tempPopulation=[Population,Offspring,AS];
                    [~,ia,~]=unique(tempPopulation.objs,'rows');
                    tempPopulation=tempPopulation(ia);
%                     childIndex=Problem.N+1:length(tempPopulation);
                      tempZMax=zMax;
                      tempZMin=zMin;
                      index=find(tempZMax-tempZMin<10.^-20);
                      if(length(index)~=0)
                            tempZMax(index)=1;
                            tempZMin(index)=0;
                      end
                     range= tempZMax-tempZMin;
                     fuzzy=5*10.^-4*(min(range,1));
                     [FrontNo,MaxFNo] = NDSortTheta(tempPopulation.objs,min(Problem.N,length(tempPopulation)),fuzzy);
%                      [FrontNo,MaxFNo] = NDSort(tempPopulation.objs,min(Problem.N,length(tempPopulation)));
%                      AddToAS=find(FrontNo(childIndex)==1);
%                      AddToAS=offspring(FrontNo(childIndex)==1);
                      FristIndex=find(FrontNo==1);
                      normOfFirst=sum((tempPopulation(FristIndex).objs-zMin).^(2),2).^(1/2);
                      deleteIndex=normOfFirst>10*mean(normOfFirst);
                      FrontNo(FristIndex( deleteIndex))=2;
                       AS=tempPopulation(FrontNo==1);
                     
                     
                      objs=AS.objs;
                      [zMax,~] = max(objs,[],1);
                      tempZMax=zMax;
                      tempZMin=zMin;
                      index=find(tempZMax-tempZMin<10.^-20);
                      if(length(index)~=0)
                            tempZMax(index)=1;
                            tempZMin(index)=0;
                      end
                      temp=(objs-tempZMin)./(tempZMax-tempZMin);
                      
                      %% 尝试同时考虑个体的最大值与主轴的夹角来作为选择极值点的标准。
%                       cosineToLine=pdist2(temp,ones(1,M),'cosine');
%                       [~,maxIndex]=max(cosineToLine);
%                       cosineToLine(maxIndex)=[];
%                       notInMaxIndex=setdiff(1:size(temp,1),maxIndex);
%                       minCosDis=min(pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine'),cosineToLine);
%                 %       minCosDis=pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine')+cosineToLine;
%                       while length(maxIndex)<M && length(notInMaxIndex)>0
%                           [~,appendIndex]=max(minCosDis);
%                           maxIndex=[maxIndex,notInMaxIndex(appendIndex)];
%                           minCosDis(appendIndex)=[];
%                           notInMaxIndex(appendIndex)=[];
%                           minCosDis=min(minCosDis,pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine'));
%                       end
                       %% 尝试同时考虑个体的最大值与主轴的夹角来作为选择极值点的标准。
                      maxTemp=max(temp,[],2);
% %                       cosineToLine=pdist2(temp,ones(1,M),'cosine');
                      cosineToLine=pdist2(temp,ones(1,M),'cosine').*maxTemp;
                      maxTemp=max(temp,[],2);
%                       % 
                      [~,maxIndex]=max(cosineToLine);
                      cosineToLine(maxIndex)=[];
                      notInMaxIndex=setdiff(1:size(temp,1),maxIndex);
                      minCosDis=min(pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine').*maxTemp(notInMaxIndex),cosineToLine);
                %       minCosDis=pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine')+cosineToLine;
                      while length(maxIndex)<M && length(notInMaxIndex)>0
                          [~,appendIndex]=max(minCosDis);
                          maxIndex=[maxIndex,notInMaxIndex(appendIndex)];
                          minCosDis(appendIndex)=[];
                          notInMaxIndex(appendIndex)=[];
                          minCosDis=min(minCosDis,pdist2(temp(notInMaxIndex,:),temp(maxIndex(end),:),'cosine').*maxTemp(notInMaxIndex));
                      end
                     
                     ExtremeVectors=temp(maxIndex,:);
                     if(size(GNGnet.NodeP,1)>M)
                         [~,EMIndex]=max(1-pdist2(ExtremeVectors,GNGnet.NodeP,'cosine'),[],2);
                         EM=GNGnet.NodeP(EMIndex,:);
                     else
%                          EM=eye(M)+10.^-6;
                        EM=ExtremeVectors;
                     end
                     EM=0.99*ExtremeVectors./(sum(ExtremeVectors,2))+0.01*(EM);
%                      
                     ChooseIndex = ArchiveUpdate_MD(AS.objs,ArchiveSize1,GNGnet,zMin,zMax,EM);
                     AS=AS(ChooseIndex);

                    
                    
%                     [AS,~] = ArchiveUpdate_MD([AS;Population.objs],length(Population),ArchiveSize1,ArchiveSize2,GNGnet,zMin,zMax);
                    nAS = length(AS);
                    
                    % GNG Update (and Algorithm 3)
%                     GNGnet.maxNode = min(ceil(1*Problem.N),floor(nAS/2)); % paramter reset 
                     GNGnet.maxNode = min(ceil(1*Problem.N));
                        GNGnet.maxHP = 2*nAS; % paramter reset
%                     if(mod(ceil(Problem.FE/Problem.maxFE),2))
                    GNGnet = GNGUpdateByDirection(AS.objs,GNGnet,Problem.M+1,zMin);
%                     Ru=(AS.objs-tempZMin)./(tempZMax);
%                         PopulationIndex=1:Problem.N;
%                         AS
%                         [AS,~] = ArchiveUpdate_MD([AS;Population.objs],length(Population),ArchiveSize1,ArchiveSize2,GNGnet,zMin,zMax);
%                     end
%                     theta = TunePBI(GNGnet,eps); % Tune theta in PBI function
%                     theta = theta+1-min(theta);
                     %% Environmental Selection
%                     theta = 1/2*TunePBI1(GNGnet,eps); % Tune theta in PBI function
%                    theta=zeros(1,size(GNGnet.NodeS,1))+0.5;
%                     [Population,FrontNo,GNGnet,zMax] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru,GNGnet,zMin,zMax,(Problem.FE/(1*Problem.maxFE)).^(3),EM);
%                      Ru1=(AS.objs-tempZMin)./(tempZMax);
%                          [Population,FrontNo,GNGnet,zMax] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru,GNGnet,zMin,zMax,0.75*(Problem.FE/(1*Problem.maxFE)).^(1),EM);
%                      [Population,FrontNo,GNGnet,zMax] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru1,GNGnet,zMin,zMax,0.75*(Problem.FE/(1*Problem.maxFE)).^(1),EM);
                       [Population,FrontNo,GNGnet,zMax] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru,GNGnet,zMin,zMax,0.75*(Problem.FE/(1*Problem.maxFE)).^(1),EM);
%                       [Population,FrontNo,GNGnet,zMax] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru,GNGnet,zMin,zMax,1*(Problem.FE/(1*Problem.maxFE)).^(1),EM);
%                         [Population,FrontNo,GNGnet,zMax] = ESelection_EM(tempPopulation,FrontNo,MaxFNo,Problem.N,Problem.M,Ru,GNGnet,zMin,zMax,theta,EM); 
%                             if(mod(ceil(Problem.FE/Problem.maxFE),10))
% %                         GNGnet = GNGUpdateByDirection(AS,GNGnet,Problem.M+1,zMin);
% %                         PopulationIndex=1:Problem.N;
% %                        
%                         [AS,~] = ArchiveUpdate_MD([AS;Population(FrontNo==1).objs],sum(FrontNo==1),ArchiveSize1,ArchiveSize2,GNGnet,zMin,zMax);
%                     end
                    
                    
                else
%                     evaTimes=0;
%                     while(evaTimes<=Problem.N)
                        %% 原始方案。
%                         MatingPool = TournamentSelection(2,2,FrontNo,crd);
%             %           对整个优化使用小生境的方法来产生子代。
%                         try
%                             Offspring  = OperatorGAhalf(Problem,Population(MatingPool));
%                         catch
%                             temp=0;
%                         end
%                             evaTimes= evaTimes+1;
%                         zMin= min([zMin;Offspring.objs],[],1);
                        %% 此后使用基于LMPFE中拟合前沿的方式来依据投影点衡量个体的多样性以及收敛性。
                        % 估计前沿，首先利用GNG网络来判断前沿的分割区域，对不同的分割区域使用对应关联的种群个体来进行拟合。
                        if(ceil(Problem.FE/Problem.N) == MaxGen - NoG+1 && trainedFlag==0)
                            PopObjs=Population.objs;
                            FirstIndex=find(FrontNo==1);
                            objs=PopObjs(FirstIndex,:);
                            [zMax,maxIndex]=max(objs,[],1);
                            
%                             AS = ArchiveUpdate_MD([AS;Population(FrontNo==1).objs],sum(FrontNo==1),ArchiveSize1,ArchiveSize1,GNGnet,zMin,zMax);
                             tempAS=AS;     
                            [Fronts,AS]=ApproximateFrontier(GNGnet,Population(FrontNo==1).objs,zMax,zMin);
%                             tempPopulation=[Population,AS];
% %                              tempPopulation=[Population,Offspring,AS];
%                             [~,ia,~]=unique(tempPopulation.objs,'rows');
%                             tempPopulation=tempPopulation(ia);
                            
%                             AS = ArchiveUpdate_MD([AS;Population(FrontNo==1).objs],sum(FrontNo==1),ArchiveSize1,ArchiveSize1,GNGnet,zMin,zMax);
%                             [tempFrontNo,~] = NDSortTheta(tempPopulation.objs,min(Problem.N,length(tempPopulation)),fuzzy);
%                             [Fronts,AS]=ApproximateFrontier(GNGnet,tempPopulation(tempFrontNo==1).objs,zMax,zMin);
%                           
%                             AS=tempPopulation(tempFrontNo==1);
%                             ChooseIndex = ArchiveUpdate_MD(AS.objs,ArchiveSize1,GNGnet,zMin,zMax,EM);
%                             AS=AS(ChooseIndex);
%                              tempAS=AS;       
%                             [Fronts,AS]=ApproximateFrontier(GNGnet,tempAS.objs,zMax,zMin);
                            %% 计算初始时的收敛度以及多样性的衡量指标。
                            [app,pp]=ProjectPoints(Fronts,Population.objs,AS,zMax,zMin);
                            [div,cov]=cpDivAndCov(Population.objs,app,zMin);
                            preDiv=mean(div);
                            preCov=mean(cov);
                            preDivs=[preDiv];
                            preCovs=[preCov];
                            theta2=0.7;
                            % theta2=abs(preCov/preDiv);
                            theta2s=[theta2];
                            trainedFlag=1;
                            
                            F1 = Population(FrontNo==1);
                            objs=F1.objs;
%                             [zMax,~]=max(objs,[],1);
                            [zMax,~] = max(AS,[],1);
                            info.zMax=zMax;
                            info.zMin=zMin;
                            info.TrueZMin=zMin;
                            
                            info.pp=pp;
                            PopObjs=(Population.objs-zMin)./(zMax-zMin);
                            info.cov=(app-sum(PopObjs,2))./(app);
                            info.dis=squareform(pdist(pp));
                            [sortDis,order]=sort(info.dis,2,'ascend');
                            info.sortDis=sortDis;
                            info.order=order;
                            info.Fno=FrontNo;
                            temp=0;
                            crd=zeros(1,Problem.N);
                        end
                      
                        %% 得到拟合前沿后，利用投影点进行筛选(但问题就在于估计不准，能与真实前沿的收敛程度达到0.1.....)，所以在由于还是使用基于支配等级优先。
                        % [Population,FrontNo,crd,EM]=ESelection_AF([Population,Offspring],Problem.N,Ruq,GNGnet,Fronts,NetLabel,theta2,zMin,EM);
                        % 版本5: 采用(miu+1)进化策略来选择。
                        zMin=info.zMin;
                        if(Problem.FE>=27073)
                            temp=0;
                        end
                        [Population,FrontNo,info,crd]=ESelection_OneByOne(Problem,crd,Population,GNGnet,Fronts,info);
                        % 版本5: 采用(miu+1)进化策略来选择。

                        [app,~]=ProjectPoints(Fronts,Population.objs,AS,zMax,zMin);
                        [div,cov]=cpDivAndCov(Population.objs,app,zMin);


                        newDiv=mean(div);
                        newCov=mean(cov);
                        [theta2,preCov,preDiv]=UpdateTheta(preCov,preDiv,newCov,newDiv);
                        preDivs=[preDivs,preDiv];
                        preCovs=[preCovs,preCov];
                        theta2s=[theta2s;theta2];
%                     end
                end
            end
            
             folder = fullfile('.\Data',class(Algorithm),class(Problem),num2str(Problem.M));
            [~,~]  = mkdir(folder);
            % folder=folder(1:end-1);
            file   = fullfile(folder,sprintf('GNGnetAndFronts_%s_%s_M%d_D%d_N%d_',class(Algorithm),class(Problem),Problem.M,Problem.D,Problem.N));
            runNo  = 1;
            while exist([file,num2str(runNo),'.mat'],'file') == 2
                runNo = runNo + 1;
            end
            nodes=GNGnet.NodeS;
            save([file,num2str(runNo),'.mat'],'GNGnet','tempAS','Fronts','preDivs','preCovs','theta2s');
            temp=0;
            
        end
    end
    
end